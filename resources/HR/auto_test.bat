For /R IITM %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -k -i %%? -o %%~dp?%%~n?_exact_exp_output.mxt
)
forfiles /S /M *.mxt /C "cmd /c rename @file @fname.txt"