def al_conductivity(T,r0=2.97e-9):
  # AL_CONDUCTIVITY computes k(T) given T and residual resistivity r0.
  #
  # k = al_conductivity(T,r0) gives aluminum conductivity k(T) in [W/m/K]
  # using the results of Woodcraft, Cryogenics 45 (6) 421-431 (2005), where
  # r0 is the residual resistivity in [Ohm m] in the limit of T = 0 K.
  #
  # Note that this does not give valid results for superconducting Al. The
  # T_c for various alloys of Al may vary from 400 mK to 1.4 K.
  #
  # k = al_conductivity(T) defaults to the conductivity values consistent
  # with measured BICEP3 Al 1100-O (11.15 W/m/K @ 1 K, or r0 = 2.197e-9).
  #
  # For 6061-T6, r0 = 1.0125e-8. A very good 1100-O has r0 upwards of 8e-10.
  #
  # Adapted from a matlab script by Ki Won Yoon.

  L0 = 2.45e-8           #g Lorenz number [W Ohm / K^2]
  r_rt_pure = 2.43e-8    #g Room temp resistivity of pure Al [Ohm m]
  beta = r0 / L0
  RRR_p = r_rt_pure / r0 #g RRR w.r.t. r_rt_pure

  P1 = 4.716e-8
  P2 = 2.446
  P3 = 623.6
  P4 = -0.16
  P5 = 130.9
  P6 = 2.5
  P7 = 0.8168

  a1 = 2.958e-8
  b1 = 0.129
  a3 = 925.4
  b3 = -0.167

  P1R = min([a1 * RRR_p^b1, P1])
  P3R = max([a3 * RRR_p^b3, P3])

  W0 = beta ./ T

  Wi = P1R.*T.^P2 ./ (1 + P1R.*P3R.*T.^(P2+P4).*exp(-1.*((P5./T).^P6)))

  Wi0 = P7 .* Wi .* W0 ./ (Wi + W0)

  k = 1 ./ (W0 + Wi + Wi0)

  return k
