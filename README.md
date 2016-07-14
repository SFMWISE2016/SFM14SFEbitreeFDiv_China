# SFM14SFEbitreeFDiv_China
Final project #3 by group #14 for SFM2016 course
"Application of the Quantlet of the Binomial tree to CSI 300 Index Option (China)

```yaml
Name of QuantLet: SFEbitreeFDiv_China

Description: 'Computes European option prices using a binomial tree for assets with fixed amount dividends.'

Keywords: 'binomial, tree, asset, call, put, option, option-price, european-option, dividends, financial, black-scholes'

See also: SFEbitreePDiv, SFEbitreeCDiv, SFSbitreeNDiv, SFEbitreeFDiv

Author: Awdesch Melzer, Piotr Majer, Ying Chen, Anait Agadzhanian, Jiang Mengjie, Wouter Biegstraaten

Submitted: Fri, July 15 2015 by Anait Agadzhanian

Input: 
- s0 : Stock Price
- k : Exercise Price
- i : Interest Rate
- sig : Volatility
- t : Time to Expiration
- n : Number of Intervals
- type : 0 is American, 1 is European
- flag : 1 is call, 0 is Put
- div: 0

Output: binomial trees and price of option

Example: 'User inputs parameters [s0, k, i, sig, t, n, type] like [2.136, 2, 0.0305, 0.13, 0.153424658, 5, 1， 0], then call price is shown.'

```

```r
# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()


s0 = 2.136  # Stock price
k = 2  # Exercise price
i = 0.0305  # Rate of interest
sig = 0.13  # Volatility
t = 0.153424658  # Time to expiration
n = 5  # Number of intervals
type = 1  # 0 is American/1 is European
flag = 1  # 1 is call/0 is put
div = 0  # Contionous dividend in percentage

# Check conditions

if (div < 0) {
    print("SFEBiTree: Dividend must be nonnegative! Please input again. div=")
    div = scan()
}

if (s0 <= 0) {
    print("SFEBiTree: Price of Underlying Asset should be positive! Please input again. s0=")
    s0 = scan()
}
if (k < 0) {
    print("SFEBiTree: Exercise price couldnot be negative! Please input again. k=")
    k = scan()
}
if (sig < 0) {
    print("SFEBiTree: Volatility should be positive! Please input again. sig=")
    sig = scan()
}
if (t <= 0) {
    print("SFEBiTree: Time to expiration should be positive! Please input again. t=")
    t = scan()
}
if (n < 1) {
    print("SFEBiTree: Number of steps should be at least equal to 1! Please input again. n=")
    n = scan()
}


if (div < 0) {
    print("SFEBiTree: Dividend must be nonnegative! Please input again. div=")
    div = scan()
}

# Main computation
dt = t/n  # Interval of step
u = exp(sig * sqrt(dt))  # Up movement parameter u
d = 1/u  # Down movement parameter d
b = i  # Costs of carry
p = 0.5 + 0.5 * (b - sig^2/2) * sqrt(dt)/sig  # Probability of up movement
tdivn = floor(tdiv/t * n - 1e-04) + 1
s = matrix(1, n + 1, n + 1) * s0
un = rep(1, n + 1) - 1
un[n + 1] = 1
dm = t(un)
um = matrix(0, 0, n + 1)
j = 1

while (j < n + 1) {
    d1 = cbind(t(rep(1, n - j) - 1), t((rep(1, j + 1) * d)^(seq(1, j + 1) - 1)))
    dm = rbind(dm, d1)  # Down movement dynamics
    u1 = cbind(t(rep(1, n - j) - 1), t((rep(1, j + 1) * u)^((seq(j, 0)))))
    um = rbind(um, u1)  # Up movement dynamics
    j = j + 1
}

um = t(rbind(un, um))
s = s[1, 1] * um * dm  # Stock price development
print("Stock_Price")
print(s)
s = s[nrow(s):1, ]  # Rearangement
opt = matrix(0, nrow(s), ncol(s))

if ((flag == 1) && (type == 1)) {
    # Option is a European call
    opt[, n + 1] = pmax(s[, n + 1] - k, 0)  # Determine option values from prices
    loopv = seq(n, 1)
    for (j in loopv) {
        l = seq(1, j)
        # Probable option values discounted back one time step
        discopt = ((1 - p) * opt[l, j + 1] + p * opt[l + 1, j + 1]) * exp(-b * dt)
        # Option value
        opt[, j] = rbind(t(t(discopt)), t(t(rep(0, n + 1 - j))))
    }
    European_Call_Price = opt[nrow(opt):1, ]
    print(European_Call_Price)
    print(" ")
    print("The price of the option at time t_0 is")
    print(European_Call_Price[n + 1, 1])
}

if ((flag == 0) && (type == 1)) {
    # Option is a European put
    opt[, n + 1] = pmax(k - s[, n + 1], 0)  # Determine option values from prices
    loopv = seq(n, 1)
    for (j in loopv) {
        l = seq(1, j)
        # Probable option values discounted back one time step
        discopt = ((1 - p) * opt[l, j + 1] + p * opt[l + 1, j + 1]) * exp(-b * dt)
        # Option value
        opt[, j] = rbind(t(t(discopt)), t(t(rep(0, n + 1 - j))))
    }
    European_Put_Price = opt[nrow(opt):1, ]
    print(European_Put_Price)
    print(" ")
    print("The price of the option at time t_0 is")
    print(European_Put_Price[n + 1, 1])
}
```
