import os
import numpy as np
from astropy.table import Table
from astroquery.gaia import Gaia


def dr2source_edr3source(dr2_source_id):
    """
    This function will find the best EDR3 source id match provided in the table `gaiaedr3.dr2_neighbourhood`

    Documentation for `gaiaedr3.dr2_neighbourhood` available here: https://gaia.aip.de/metadata/gaiaedr3/dr2_neighbourhood/

    :return: row-matched DR3 source id, otherwise 0
    """
    neighbourhood_table = "gaiaedr3.dr2_neighbourhood"
    t = Table({'source_id': np.unique(dr2_source_id)})

    # try to login if any
    if os.stat("gaia_credential").st_size != 0:
        Gaia.login(credentials_file='gaia_credential')

    # launching job at https://gea.esac.esa.int/archive/
    job = Gaia.launch_job_async(
        f"""
        select n.* 
        from {neighbourhood_table} as n 
        inner join tap_upload.my_table as m on m.source_id = n.dr2_source_id
        """,
        upload_resource=t,
        upload_table_name="my_table")

    result = job.results

    dr2source_id = np.array(result["dr2_source_id"])
    dr3source_id = np.array(result["dr3_source_id"])
    angular_distance = np.array(result["angular_distance"])
    magnitude_difference = np.array(result["magnitude_difference"])
    idx = (dr2source_id != dr3source_id)  # indices where those source id changed in dr3
    idx_inverse = np.where(~idx)[0]  # indices where those source id unchanged in dr3
    idx = np.where(idx)[0]

    # this loop makes sure if there is an exact source_id suggestion in DR3, use that and delete all other suggestions
    while len(np.intersect1d(dr2source_id[idx], dr2source_id[idx_inverse])) != 0:  # make sure its empty
        intersec, idx1, idx2 = np.intersect1d(dr2source_id[idx], dr2source_id[idx_inverse], return_indices=True)

        dr2source_id = np.delete(dr2source_id, idx[idx1])
        dr3source_id = np.delete(dr3source_id, idx[idx1])
        angular_distance = np.delete(angular_distance, idx[idx1])
        magnitude_difference = np.delete(magnitude_difference, idx[idx1])

        idx = (dr2source_id != dr3source_id)
        idx_inverse = np.where(~idx)[0]  # indices where those source id unchanged in dr3
        idx = np.where(idx)[0]

    # deal with uniqueness of the source id (cases where multiple DR3 id matched to one DR2 id)
    unique, unique_idx, unique_counts = np.unique(dr2source_id, return_index=True, return_counts=True)
    unique_m1_idx = np.where(unique_counts > 1)  # those cases with multiple DR3 id
    bad_dr2id_matches = []  # record bad dr2 id
    bad_dr2id_idx = np.array([], dtype=int)  # bad dr2 id indices

    for _source_id in unique[unique_m1_idx]:
        _idx = np.where(dr2source_id == _source_id)[0]
        if np.argmin(angular_distance[_idx]) == np.argmin(np.abs(magnitude_difference[_idx])):
            _idx = np.delete(_idx, np.argmin(angular_distance[_idx]))
            bad_dr2id_matches.append(dr2source_id[_idx])
            bad_dr2id_idx = np.concatenate([bad_dr2id_idx, _idx])
        else:  # assume no good match, put the whole thing to bad matches array
            bad_dr2id_matches.append(dr2source_id[_idx])
            bad_dr2id_idx = np.concatenate([bad_dr2id_idx, _idx])

    dr2source_id = np.delete(dr2source_id, bad_dr2id_idx)
    dr3source_id = np.delete(dr3source_id, bad_dr2id_idx)

    # perform row-matching
    index = np.argsort(dr2source_id)
    sorted_x = dr2source_id[index]
    sorted_index = np.searchsorted(sorted_x, dr2_source_id)
    yindex = np.take(index, sorted_index, mode="clip")
    mask = dr2source_id[yindex] != dr2_source_id
    result = np.ma.array(yindex, mask=mask)
    matched = dr3source_id[result.filled(0)]
    matched[mask] = 0

    return matched
