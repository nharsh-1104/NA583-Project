SUDDEN ROAD CHANGE - PESC ABS (ORIGINAL CONTROLLER LEFT UNTOUCHED)
==================================================================

Run this one MATLAB script:
    run_pesc_road_change

FILES IN THIS PACKAGE
  ORIGINAL, BYTE-FOR-BYTE UNMODIFIED:
    abs_plant.m
    pesc.m
    slip_controller.m
    run_pesc_constant_mu.m

  NEW FILES:
    abs_plant_roadchange.m    Variable road friction: dry -> ice at 45 m
    pesc_roadchange.m         Road-change detector + reset/freeze actions
    run_pesc_road_change.m    Two-controller comparison / plots
    README_sudden_road_change.txt

VERY IMPORTANT: ORIGINAL CONTROLLER SETTINGS ARE NEVER RETUNED
  Original p.es.a     = 0.02
  Original p.es.w     = 10
  Original p.es.wh    = 2
  Original p.es.wl    = 3
  Original p.es.k     = 8 (NOT 2)
  Original p.es.phi   = 0
  Original initial lambda_hat = 0.10
  Original p.c        = 20

  NO startup adaptation hold is applied to either controller.
  NO additional reference saturation or dither clipping is applied.
  The original PESC function pesc.m is directly called in BOTH cases.
  The modified controller only adds friction-drop detection and temporary
  recovery (reset slip estimate/filters and freeze adaptation for 0.30 s).

ROAD SCENARIO (same for BOTH controllers)
  Dry road until distance 45 m: peak mu=0.80, peak slip=0.12
  Icy road thereafter         : peak mu=0.16, peak slip=0.05
  Both conditions are illustrative analytical friction curves.

SIMULATION
  2-ms sampled PESC (explicit Euler) and RK4 wheel/vehicle plant.
  Same integration algorithm, initial conditions, and road for both.
  120 km/h initial speed, ends at 1 m/s or 25 s.
  The main script saves PNG comparison and a MAT data file.

CONTROLLER CHANGES (ONLY in pesc_roadchange.m)
  Detect a measured mu drop >0.08 in one 2-ms sample while change in
  measured slip is <=0.015. The detector is enabled at 0.5 s; this does
  NOT freeze the extremum-seeking integrator before detection.
  On a detected drop:
     lambda_hat = min(lambda_hat, 0.075) (conservative fallback)
     high-pass LPF state = current measured friction
     demodulated LPF state = 0
     pause gradient adaptation for 0.30 s, then resume original PESC
     suppress another detection for 0.80 s

CAVEATS
  * Original controller file is unmodified, including its initial
    transient and current gain; that may cause poor initial dry-road slip.
  * Measured friction is the exact plant signal, not sensor-estimated mu.
    The detector is a proof-of-concept, not yet robust to road sensor noise.
  * The original constant-road demo uses ode45; this fair comparison uses
    a common sampled implementation of its exact PESC derivative.
  * The detector uses NO true road-change location or actual ice slip.
  * MATLAB was not available in the creation environment; numerical output
    should be verified by running the script in MATLAB.
  * All statements about relative stopping distance require verification
    from the actual run, rather than assuming the modified controller wins.
