import pandas as pd

# Read the CSV file
fincham_df = pd.read_csv('Data/Fincham1992_energy.csv')
fincham_df = fincham_df.groupby(['dt']).agg({'error': ['mean', 'std']}).reset_index()
fincham_df.columns = ['Time_Fincham', 'Fincham', 'Fincham_error']


omelyan_df = pd.read_csv('Data/Omelyan1998_energy.csv')
omelyan_df = omelyan_df.groupby(['dt']).agg({'error': ['mean', 'std']}).reset_index()
omelyan_df.columns = ['Time_Omelyan', 'Omelyan', 'Omelyan_error']


spiral_df = pd.read_csv('Data/delValle2023_energy.csv')
spiral_df = spiral_df.groupby(['dt']).agg({'error': ['mean', 'std']}).reset_index()
spiral_df.columns = ['Time_SPIRAL', 'SPIRAL', 'SPIRAL_error']


merged_df = pd.merge(fincham_df, omelyan_df, how='outer', left_on='Time_Fincham', right_on='Time_Omelyan', suffixes=('_Fincham', '_Omelyan'))
merged_df = pd.merge(merged_df, spiral_df, how='outer', left_on='Time_Fincham', right_on='Time_SPIRAL', suffixes=('_Merged', '_SPIRAL'))
merged_df.to_csv('Data/Energy.csv', index=False)
