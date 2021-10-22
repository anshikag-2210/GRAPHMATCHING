break>.\smfq_stats_master_10_5.csv
break>.\smfq_stats_master_20_5.csv
break>.\smfq_stats_master_100_5.csv
break>.\smfq_stats_master_1000_5.csv
break>.\smfq_stats_random_10_5.csv
break>.\smfq_stats_random_20_5.csv
break>.\smfq_stats_random_100_5.csv
break>.\smfq_stats_random_1000_5.csv
break>.\smfq_stats_shuffle_10_5.csv
break>.\smfq_stats_shuffle_20_5.csv
break>.\smfq_stats_shuffle_100_5.csv
break>.\smfq_stats_shuffle_1000_5.csv

For /R HR\master\n1_1000_n2_10_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_master_10_5.csv
)
For /R HR\master\n1_1000_n2_20_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_master_20_5.csv
)
For /R HR\master\n1_1000_n2_100_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_master_100_5.csv
)
For /R HR\master\n1_1000_n2_1000_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_master_1000_5.csv
)

For /R HR\random\n1_1000_n2_10_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_random_10_5.csv
)
For /R HR\random\n1_1000_n2_20_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_random_20_5.csv
)
For /R HR\random\n1_1000_n2_100_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_random_100_5.csv
)
For /R HR\random\n1_1000_n2_1000_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_random_1000_5.csv
)

For /R HR\shuffle\n1_1000_n2_10_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_shuffle_10_5.csv
)
For /R HR\shuffle\n1_1000_n2_20_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_shuffle_20_5.csv
)
For /R HR\shuffle\n1_1000_n2_100_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_shuffle_100_5.csv
)
For /R HR\shuffle\n1_1000_n2_1000_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o .\test_output.txt >> .\smfq_stats_shuffle_1000_5.csv
)


