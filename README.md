# SOFtools

A small toolbox aimed to design static output feedback with an Hinfinity performance objective. 

The synthesis problem involves solving a nonconvex optimization problem. The toolbox implements several approximation methods proposed in:
1. M. Mattei, ‘Robust multivariable PID control for linear parameter varying systems’, Automatica, vol. 37, no. 12, pp. 1997–2003, 2001
2. C.Crusius and A. Trofino, ‘Sufficient LMI conditions for output feedback control problems’, IEEE Trans. on Automatic Control, vol. 44, no. 5, pp. 1053–1057, 1999.
3. E. Prempain, “Static output feedback stabilisation with Hinf performance for a class of plants,” Syst. Control Lett., vol. 43, no. 3, pp. 159–166, Jul. 2001.
4. L. El Ghaoui, F. Oustry, and M. AitRami, ‘A cone complementarity linearization algorithm for static output-feedback and related problems’, IEEE Trans. on Automatic Control, vol. 42, no. 8, p. 1171, 1997.

*Use:*
`[Ksof,gopt] = sofsyn(sys,ios, [method], [opt])`
argument within [] are optional
*Inputs:*
- `sys`: augmented plant (LTI or LPV)
- `ios`: set the measured outputs and the control inputs
- `method`: mattei (default)/crucius/prempain/lowrank/systune
*Outputs:*
- `Ky`: static output gain
- `gopt`: optimal performance level

Examples of how to use it can be found in the folder **demo**.


