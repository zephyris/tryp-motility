FOR /L %%A IN (1,1,30) DO (
  ECHO Starting thread %%A
  START C:\ImageJ\ParticleTrackingv4.bat %%A
  TIMEOUT 1
)