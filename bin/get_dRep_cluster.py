#!/usr/bin/env python

import pandas as pd
import sys



with open(sys.argv[3],'w') as outF:
#def read_data_tables(wdb,cdb):
    cluster_info = {}
    #wdb = data_tables_dir / "Wdb.csv"
    #cdb = data_tables_dir / "Cdb.csv"
    wdb_df = pd.read_csv(sys.argv[1], index_col=None, header=0)
    cdb_df = pd.read_csv(sys.argv[2], index_col=None, header=0)
    for _, row in wdb_df.iterrows():
        cluster_info.setdefault(row["genome"], set())
        cluster_id = row["cluster"]
        cluster_info[row["genome"]].update(
            cdb_df.loc[
                (cdb_df["secondary_cluster"] == cluster_id)
                & (cdb_df["genome"] != row["genome"]),
                "genome",
            ].tolist()
        )

    for a in cluster_info:
        outF.write(str(a)+'\t'+';'.join(cluster_info[a])+'\n')


