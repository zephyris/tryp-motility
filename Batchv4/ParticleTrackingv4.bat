for %%v in (*.czi) do C:\ImageJ\jre\bin\java.exe -Xmx16g -jar C:\ImageJ\ij.jar -ijpath C:\ImageJ -batch C:\ImageJ\ParticleTrackingv4.ijm "%%v*%cd%"
EXIT