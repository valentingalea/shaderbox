set NAME=hlsltoy
set VER=1.0.0.0
set PUB="Valentin Galea"
set URL=\\VALENTIN-PC\setup\%VER%\%NAME%.application"
REM "http://valentingalea.github.io/shaderbox/"

mage -New Application -Processor x86 -ToFile %VER%\%NAME%.exe.manifest -name "%NAME%" -Version %VER% -FromDirectory %VER%
mage -New Deployment -Processor x86  -ToFile %VER%\%NAME%.application -Install true -Publisher %PUB% -ProviderUrl %URL% -AppManifest %VER%\%NAME%.exe.manifest -AppCodeBase %VER%\%NAME%.exe.manifest
