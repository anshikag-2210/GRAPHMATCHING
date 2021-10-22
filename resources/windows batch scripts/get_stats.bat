break>.\popular1_stats_master_10_5.csv
break>.\popular1_stats_master_20_5.csv
break>.\popular1_stats_master_100_5.csv
break>.\popular1_stats_shuffle_10_5.csv
break>.\popular1_stats_shuffle_20_5.csv
break>.\popular1_stats_shuffle_100_5.csv
break>.\popular1_stats_random_20_5.csv
break>.\popular1_stats_random_100_5.csv
For /R HRLQ\master\n1_1000_n2_10_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_master_10_5.csv
)
For /R HRLQ\master\n1_1000_n2_20_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_master_20_5.csv
)
For /R HRLQ\master\n1_1000_n2_100_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_master_100_5.csv
)
For /R HRLQ\shuffle\n1_1000_n2_10_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_shuffle_10_5.csv
)
For /R HRLQ\shuffle\n1_1000_n2_20_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_shuffle_20_5.csv
)
For /R HRLQ\shuffle\n1_1000_n2_100_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_shuffle_100_5.csv
)
For /R HRLQ\random\n1_1000_n2_20_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_random_20_5.csv
)
For /R HRLQ\random\n1_1000_n2_100_k_5 %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -z -i %%? -o .\test_output.txt >> .\popular1_stats_random_100_5.csv
)
