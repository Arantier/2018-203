@echo off
set name=scherbakovdv
cls
echo ���� ���������� �� �ࠪ⨪� �������!
set /p num="������ ����� ���� ��� �஢�ન, 0 ��� �������樨, 10 ��� �ਢ������ ��४�ਨ � 祫����᪨� ���:"
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