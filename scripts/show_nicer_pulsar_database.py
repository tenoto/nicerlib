#!/usr/bin/env python

import os 
import pandas as pd

csvfile = '%s/scripts/data/nicer_pulsar_database_te180406.csv' % os.getenv('NICER_SOFT_PATH')
df = pd.DataFrame.from_csv(csvfile)
print(df)

