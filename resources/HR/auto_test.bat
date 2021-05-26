For /R IITM %%? IN (*.txt) do (
C:\Users\SAI\source\repos\GraphMatching\build\graphmatching -B -l -i %%? -o %%~dp?%%~n?_new_output.mxt
)
forfiles /S /M *.mxt /C "cmd /c rename @file @fname.txt"