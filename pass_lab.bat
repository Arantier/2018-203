@echo off
set name=scherbakovdv
cls
echo ���� ���������� �� �ࠪ⨪� �������!
set /p num="������ ����� ���� ��� �஢�ન, 0 ��� �������樨, 10 ��� �ਢ������ ��४�ਨ � 祫����᪨� ���:"
if %num%==0 (	
	"C:\Program Files\Git\mingw64\bin\make"
	del *.o
) else if %num% GEQ 1 if %num% LEQ 9 (
vvm %name% %num%
) else if %num%==10 (	del vvm.exe
) else echo You obosralsya
pause