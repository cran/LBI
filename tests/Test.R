library(LBI)
# Extreme value test

LIbin(0, 1)
LIbin(0, 2)
LIbin(0, 10)
LIbin(0, 1e308)
LIbin(1, 1)
LIbin(10, 10)
LIbin(1e308, 1e308)

LIpois(0)
LIpois(1)
LIpois(2)
LIpois(3)

#LInorm
LInorm(c(0, 1, 2))

RDLI(0, 1, 0, 1)
RDLI(0, 1, 1, 1)
RDLI(0, 10, 0, 10)
RDLI(0, 10, 10, 10)
RDLI(1, 1, 0, 1)
RDLI(1, 1, 1, 1)
RDLI(10, 10, 0, 10)
RDLI(10, 10, 10, 10)
