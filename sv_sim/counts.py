import os as os
import pandas as pd
import pyranges as pr
from glob import glob

result = [y for x in os.walk("/Data/Analyses/2022/202205_SV-SIM/") for y in glob(os.path.join(x[0], '*summary.bed'))]

svs = pd.DataFrame(columns=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'],)

for fi in result:
    df = pr.read_bed(fi, as_df=True)
    svs = svs.append(df)

# svs.set_axis({'Chromosome', 'Start', 'End', 'Name', 'Insert', 'Rand'}, axis=1, inplace=True)

svs['Name'].value_counts()
