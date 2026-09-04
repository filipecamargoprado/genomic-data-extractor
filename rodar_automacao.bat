@echo off
title Automacao MMseqs2
color 0A

echo =======================================================
echo        INICIANDO AUTOMACAO - AGRUPAMENTO MMSEQS2
echo =======================================================
echo.

cd /d "%~dp0"
python agrupar_gz_ncbi.py

echo.
echo =======================================================
echo        PROCESSO FINALIZADO!
echo =======================================================
pause
