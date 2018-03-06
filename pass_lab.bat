@echo off
set name=scherbakovdv
cls
echo Добро пожаловать на практику Жалнина!
set /p num="Введите номер лабы для проверки, 0 для компиляции, 10 для приведения директории в человеческий вид:"
if %num%==0 (	
	REM cd .\2018-203\
	set way="C:\Program Files\Git\mingw64\bin\"
	if not exist %way%make (
		set way="C:\Program Files (x86)\CodeBlocks\MinGW\bin\"
		%way%mingw32-make
	) else %way%make
	del *.o
	REM move vvm.exe ..\
) else if %num% GEQ 1 if %num% LEQ 9 (
vvm %name% %num%
) else if %num%==10 (	del vvm.exe
) else echo You obosralsya
pause