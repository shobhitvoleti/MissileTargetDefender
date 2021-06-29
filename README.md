# MissileTargetDefender
Missile Target Defender Problem
I have looked a solution to the Missile Target Defender game: A spin on the 2 player pursuit evasion game, this 3 player game includes an additional player i.e the defender that tries to stop the pursuer from intercepting the target. On inspecting the literature, most methods split the timeline into 2 parts. The first is a game between the missile, defender and the target and the second part is a game solely between the pursuer and evader. This is effected by considering two different cost functions or considering an impulse function. Based on the terminal weights, game has been classified as between decisive players, bad defender, bad pursuer and solution to each case has been obtained analytically. 


All of the solutions I've looked at considers the reduced order problem through terminal projection. The problem has also been solved for the case when the controls of each player is constrained. Sufficiency conditions for the existence of a solution have been derived. This problem has been solved analytically by Rusnak et al. as the lady, the bandit and the body guard game in literature. 

However, I've adopted the solution obtained from the Differential Riccati equation of the modified problem. As in literature, the optimisation problem has been modified and expressed using a state vector constructed from projected states derived through the state transition matrices. The control for each player has been obtained by integrating backwards in time based on terminal boundary conditions twice for two sepaprate time periods. A sufficient condition to ensure saddle point solution to exist has been stated and the system has been simulated on MATLAB. The defender is assumed to stop acting after the first duration of time. Assumption has been made that the missile, target and defender are decisive and their control actions are unbounded and simulation has been designed as such. 


The simulations have been performed on MATLAB and the .m file has been provided along with a detailed report
