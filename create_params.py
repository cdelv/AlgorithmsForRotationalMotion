# python3 create_params.py > Chess_Batch_Params.txt

print('seed', 'dt', 'integration_type')
for algorithm in ["Fincham1992", "Omelyan1998", "delValle2023"]:
	for dt in sorted([1.1e-07, 1.93117038e-07, 3.06069879e-07, 4.96369887e-07, 8.96369887e-07, 1.21848613e-06, 2.89675556e-06, 4.85088067e-6, 8.85088067e-06, 1.93117038e-05]): # 7.68812776e-8
		for seed in range(50):
			print(seed, dt, '"'+algorithm+'"')
