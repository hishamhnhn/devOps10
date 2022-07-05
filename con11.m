mo = 1.5
d = 2.2 
kc = 20.5 
F = 1.5
G = tf([1.5],[1.5 2.2 20.5])
rlocus(G)
step(G)


kp = 290 
kd = 15
Ki = 1
cont = pid(kp,ki,kd)

ref1 = feedback(cont*G,1)

t = 0:0.01:2.5;



step(ref1,t)
