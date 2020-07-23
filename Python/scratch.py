#%%
Reff = 3; a = 135; p = 0.5; num_pts = 4000

a   = np.radians(a)
asp = p * a / 2
Rp  = 1 / 2 / np.sqrt(asp)
sp  = np.sqrt(2 * asp)

s0  = 2 * sp + Rp * a * (1 - p)
scale = a / (2 * sp * (s0 - sp))
if p == 0: 
    s0 = a
    scale = a / s0

s = np.linspace(0, s0, num_pts)
K = np.zeros(num_pts)
if p == 0: K += 1
else:
    for i in range(len(K)):
        if   s[i] <= sp         : K[i] = 2 * s[i]
        elif sp < s[i] < s0 - sp: K[i] = 2 * sp
        elif s0 - sp < s[i] < s0: K[i] = 2 * (s0 - s[i])
K *= scale / Reff
s *= Reff

ds = s[1] - s[0]
φ = cumtrapz(K * ds)
x = np.cumsum(ds * np.cos(φ))
y = np.cumsum(ds * np.sin(φ))
x = np.concatenate([[0], x])
y = np.concatenate([[0], y])

# middle = int((num_pts - 1) / 2)
# η = Reff / (y[middle] + x[middle] / np.tan(a / 2))

# Rp *= Reff * η
# sp *= Reff * η
# s0 = 2 * sp + Rp * a * (1 - p)
# scale = a / (2 * sp * (s0 - sp))
# s = np.linspace(0, s0, num_pts)
# K = np.zeros(num_pts)
# if p == 0: K += Reff
# else:
#     for i in range(len(K)):
#         if   s[i] <= sp         : K[i] = 2 * s[i]
#         elif sp < s[i] < s0 - sp: K[i] = 2 * sp
#         elif s0 - sp < s[i] < s0: K[i] = 2 * (s0 - s[i])
# K *= scale
# ds = s[1] - s[0]
# φ = cumtrapz(K * ds)
# x = np.cumsum(ds * np.cos(φ))
# y = np.cumsum(ds * np.sin(φ))
# x = np.concatenate([[0], x])
# y = np.concatenate([[0], y])

pp(np.array([x,y]).T)
plt.show()
s3, K3 = curvature(x, y)
plt.plot(s, K)
plt.plot(s3, K3)
plt.show()

# %%
