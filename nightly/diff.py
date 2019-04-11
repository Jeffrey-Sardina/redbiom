def process(study):
    id_str, artifacts = study.split(',', 1)

    #Id processing
    id_num = int(id_str.replace('(', ''))

    #Artifact processing
    artifacts_set = set(artifacts.replace('[', '').replace(']','').replace(' ','').split('),'))
    artifacts_data = tuple(tuple(a.strip().replace('(', '').replace(')', '').replace("u'", "").replace("'", "").split(',')) for a in artifacts_set)
    artifacts_data = tuple(tuple(a for a in artifacts_data))

    return id_num, artifacts_data


def load_state_file(filename):
    try:
        with open(filename, 'r') as data:
            return [process(study) for study in data]
    except IOError:
        raise


def write_data_to_verify_integrity(filename, study_data):
    try:
        with open(filename, 'w') as out:
            for item in study_data:
                print(item, file=out)
    except IOError:
        raise


def out(data, filename, suffix='.out'):
    try:
        with open(filename + suffix, 'w') as out_file:
            print(data, file=out_file)
    except IOError:
        raise


def diff(study_data_old, study_data_new):
    '''
    Takes two data collections (as tuples), the old and the new, that represent study data from Qiita, and diffs these to see what has been
    added, deleted, and modified.

    Parameters:
        study_data_old: The old data
        study_data_new: The new data

    Returns
        The set of changes in the new file relative to the old file, as a tuple of sets: (added, deleted, modified).
        The added and deleted sets contain the study ids of what was added or deleted; the modified set contains the study id, artifact id, 
        and artifact data that was changed.
    '''
    ids_old = {sd[0] for sd in study_data_old}
    ids_new = {sd[0] for sd in study_data_new}

    '''
    We start by processing study IDs
        If and only if it has been added, it is in the new set but now the old one
        If and only if it had been deleted, it is in the old set but not the new one
        If it has been modified, it must be in both sets. Items in both sets are not necessarily modified, only potentially.
    I used a single sided difference (- rather than ^) since doing ^ does not differentiate added and deleted
    '''
    added = ids_new - ids_old
    deleted = ids_old - ids_new
    potentially_updated = ids_new - added

    '''
    For those studies that may have been modified, we need to look at their data.
    Anything that has been modified will have the property that it is different in the new set relative to its entry in the old set.
        Changes at or below the artifact level are considered study modifications.
    We do not care where the changes are in the old set, but rather only where they are in the new set. So we can take a sided difference.
    I did not turn study_data_old and study_data_new into sets and take their difference since that would make it harder to identify which artifact had been changed.
    '''
    ids_tags_data_old = {(sd[0], arts[0], arts[1:])
                    for sd in study_data_old
                     if sd[0] in potentially_updated
                      for arts in sd[1]}
    ids_tags_data_new = {(sd[0], arts[0], arts[1:])
                    for sd in study_data_new
                     if sd[0] in potentially_updated
                      for arts in sd[1]}
    modified = ids_tags_data_new - ids_tags_data_old

    return added, deleted, modified


if __name__ == "__main__":
    #Get the study data from files
    study_data_old = load_state_file('truncated.json')
    study_data_new = load_state_file('truncated.1.json')

    #Write out so it can be seen that the new data structures have all the same data as the old
    write_data_to_verify_integrity('reconsituted.0.json', study_data_old)
    write_data_to_verify_integrity('reconsituted.1.json', study_data_new)

    #Get what has been added, modified, and deleted and print as proof-of-concept
    added, deleted, modified  = diff(study_data_old, study_data_new)
    print('mod')
    print(modified)
    print('del')
    print(deleted)
    print('add')
    print(added)