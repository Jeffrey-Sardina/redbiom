import pandas as pd
import qiita_db
from qiita_db.util import get_artifacts_information
import biom
import joblib
import os
from collections import defaultdict
import redbiom.admin
import json
import hashlib
import random
import time
from nightly import diff


# I have adapted this code from an original file used by the Knight Lab to perform server updates.


def usable_artifacts(artifacts):
    """We only retain artifacts which are processed and are BIOM"""
    arts = []
    for art in artifacts:
        if art.visibility != 'public':
            continue
        if art.processing_parameters is None:
            continue
        if art.artifact_type != 'BIOM':
            continue
        arts.append(art)
    return arts


def get_name_and_description(art):
    """Determine the context the artifact is associated with

    The name is based off the processing parameters command name. This name
    is tagged with a stable hash based on the processing parameter values
    so that a common type of processing (e.g., pick closed reference OTUs) can
    be safely partitioned by the specific runtime parameters. The runtime
    parameters themselves are returned as a description of the context.
    """
    art_info = qiita_db.util.get_artifacts_information([art.id])[0]
    pt = art.prep_templates[0].to_dataframe()

    necessary = {'platform', 'target_subfragment', 'target_gene'}
    if not set(art_info).issuperset(necessary):
        return (None, None)

    if art_info['deprecated']:
        return (None, None)

    gene = art_info['target_gene']
    fragment = art_info['target_subfragment']
    platform = art_info['platform']

    if len(fragment) > 1:
        return (None, None)

    fragment = fragment[0]

    if fragment.lower()[0] != 'v':
        return (None, None)

    if not gene.lower().startswith('16s'):
        return (None, None)

    pp = art_info['parameters']
    alg = art_info['algorithm']
    basename = alg.split(' (')[0].replace(' ', '_')
    context_components = [basename]

    if alg.startswith(u'Pick closed-reference OTUs'):
        if 'silva' in pp[u'reference-tax'].lower():
            ref = 'SILVA'
        else:
            ref = 'Greengenes'
        context_components.append(ref)
    elif alg.startswith(u'Deblur'):
        if platform != 'Illumina':
            return (None, None)
    else:
        return (None, None)
    # Closed-Reference-Greengenes-16S-V4-100nt-a243a1'
    context_components.append(platform)
    context_components.append(gene.split()[0].upper())
    context_components.append(fragment.upper())

    parent = art.parents[0]
    if parent.processing_parameters.command.name == 'Trimming':
        trim = parent.processing_parameters.values['length']
        context_components.append("%dnt" % int(trim))

    context_components.append(art_info['algorithm_az'][:6])

    context_name = '-'.join(context_components)

    return context_name, art_info['algorithm']


def get_biom_path(art):
    """Get the BIOM path to load"""
    biompath = None
    for _, fp, fptype in art.filepaths:
        if fptype == 'biom':
            if 'deblur' in art.processing_parameters.command.name.lower():
                if fp.endswith('reference-hit.biom'):
                    biompath = fp
            else:
                biompath = fp
    return biompath


def load_study_metadata(study_id, attempt_load):
    study = qiita_db.study.Study(study_id)
    artifacts = usable_artifacts(study.artifacts())

    # if we don't have artifacts, then there isn't a reason to go on.
    if not artifacts:
        return []

    ids_tags_contexts_paths = []
    for art in artifacts:
        context, description = get_name_and_description(art)

        if context is None:
            continue

        arttag = str(art.id)

        biompath = get_biom_path(art)

        if biompath is None:
            continue

        redbiom.admin.create_context(context, description)
        ids_tags_contexts_paths.append((study_id, arttag, context, biompath))

    if attempt_load and ids_tags_contexts_paths:
        print("operating on %d" % study_id)
        df = study.sample_template.to_dataframe()
        redbiom.admin.load_sample_metadata(df)
        redbiom.admin.load_sample_metadata_full_search(df)

        preps = []
        for tag, _, _ in ids_tags_contexts_paths:
            art = qiita_db.artifact.Artifact(tag)
            pt = art.prep_templates[0].to_dataframe()

            # indexing for text search over prep is very expensive. each
            # artifact is approximately 75 seconds irrespective of the number
            # of samples. this overhead results from the the time needed to
            # stem words using nltk. SO: let's circumvent the tagging system
            # and tag outside of redbiom.

            pt.index = ['%s_%s' % (tag, i) for i in pt.index]
            preps.append(pt)

        if preps:
            prep = pd.concat(preps)
            prep.where((pd.notnull(prep)), None, inplace=True)  # replace the Nans
            redbiom.admin.load_sample_metadata(prep)
            redbiom.admin.load_sample_metadata_full_search(prep)

    # Convert to tuple since that allows the entire structure to be used
    # with set theory. This is needed since this ids_tags_contexts_paths 
    # will become an element in another list / tuple, and set operations 
    # do not work with nested lists, but do with nested tuples
    return tuple(ids_tags_contexts_paths)


def load_sample_data(tag, context, path):
    import traceback
    if not os.path.exists(path):
        print("Unable to find: %s" % path)
        return 0

    table = biom.load_table(path)
    try:
        nsdat = redbiom.admin.load_sample_data(table, context, tag)
    except ValueError:
        print("unable to load: %s, %s, %s" % (str(tag), str(context), str(path)))
        nsdat = 0
    except Exception as e:
        # there are some studies in which there are samples in the biom table
        # which lack metadata

        print(tag, context, path)
        traceback.print_exc()
        raise
    return nsdat


def delete_sample_data(study_id, tag, context, path):
    import traceback
    if not os.path.exists(path):
        print("Unable to find: %s" % path)
        return 0

    table = biom.load_table(path)
    try:
        ndeleted = redbiom.admin.delete_studies_by_id(context, tag)
    except ValueError:
        print("unable to load: %s, %s, %s" % (str(tag), str(context), str(path)))
        nsdat = 0
    except Exception as e:
        # there are some studies in which there are samples in the biom table
        # which lack metadata

        print(tag, context, path)
        traceback.print_exc()
        raise
    return ndeleted


def load_study_metadata_of_list(studies_sublist, attempt_load=True):
    """
    For a given lsit of studies, calls load_study_metadata on all
    of the studies, and returns the list that contains all of the
    results from these calls.
    This function is intended to allow a the studies list to be 
    fragmented and loaded in parallel.

    Parameters:
        studies_sublist: the sublist of the list of studies to
        process
        attempt_load: whether the study metadata should be
        attempted to be loaded into redbiom. Defaults to True
    """

    ids_tags_contexts_paths = par(joblib.delayed(load_study_metadata)(s.id, attempt_load)
                                for s in studies_sublist
                                if s.sample_template is not None)

    return ids_tags_contexts_paths


def equal_slices(super_list, sublist_number):
    """
    Given a master list, returns a list of sub-lists that all have
    length slice_length (except the last one, which will be larger
    if len(list) is not divisible by slice_length)

    Parameters:
        super_list: The super list from which sub-lists are to
        be created
        sublist_number: The number of sublists to create
    """

    slice_location_multiples = int(len(super_list) / (sublist_number))
    slices = []
    for i in range(0, sublist_number):
        start = i * slice_location_multiples
        end = (i + 1) * slice_location_multiples
        slices.append(super_list[start: end])
    slices.append(super_list[end:])

    return slices


#if __name__ == '__main__':
def update():
    # Take redbiom out of read-only mode so that it can be updated
    redbiom.admin.ScriptManager.load_scripts(read_only=False)

    # Get studies and shuffle them randomly so study sizes are
    # homogenusously distributed (on average)
    studies = list(qiita_db.study.Study.get_by_status('public'))
    random.shuffle(studies)

    # Split the list into 8 equally-lengthed parts, and load
    # those into a master list
    sub_lists = equal_slices(studies, 8)

    # Process the sub-arrays in parallel on eight CPUs in parallel
    ids_tags_contexts_paths = []
    with joblib.parallel.Parallel(n_jobs=8, verbose=50) as par:
        ids_tags_contexts_paths.append(par(joblib.delayed(load_study_metadata_of_list)(sub_lists[i])
                                    for i in range(len(sub_lists))))

    # Diff the two files so I can know what to add, remove,
    # and modify
    study_data_old_name = "study_data_old.out"
    study_data_old = diff.load_state_file(study_data_old_name)
    study_data_new = tuple(ids_tags_contexts_paths)
    added, deleted, modified = diff.diff(study_data_old, study_data_new)
    modified_study_ids = [item[0] for item in modified]

    # Now load the data into redbiom
    with joblib.parallel.Parallel(n_jobs=8, verbose=50) as par:
        nsamp = par(joblib.delayed(load_sample_data)(t, c, p)
                    for row in ids_tags_contexts_paths
                     for i, t, c, p in row
                      if i[0] in modified_study_ids or i in added)
    #And delete anything that is not in the new set
    with joblib.parallel.Parallel(n_jobs=8, verbose=50) as par:
        nsamp = par(joblib.delayed(delete_sample_data)(i, t, c, p)
                    for row in ids_tags_contexts_paths
                     for i, t, c, p in row
                      if i in deleted)

    # Put redbiom back into read-only mode for security
    redbiom.admin.ScriptManager.load_scripts(read_only=True)
    
    #Now that data has bee uploaded, write the most recent data down as the "old" data
    with open(study_data_old_name, 'w') as last_study_data:
        print(study_data_new, file=last_study_data)
    
    #Old method
    '''
    # Now load the data into redbiom
    # I need a way to delete the studies that were deleted. 
    with joblib.parallel.Parallel(n_jobs=8, verbose=50) as par:
        nsamp = par(joblib.delayed(load_sample_data)(t, c, p)
                    for row in ids_tags_contexts_paths
                     for i, t, c, p in row
                      if i in modified_study_ids or i in added)

    # Put redbiom back into read-only mode for security
    redbiom.admin.ScriptManager.load_scripts(read_only=True)
    '''
