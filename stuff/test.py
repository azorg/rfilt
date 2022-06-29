#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from math import *

r = 1.
t = 10.

v = 1. + r * t * t * 0.5
a = r * t

a = -a

if v > 0.:
  r = a * a / (2. * v)
  pass
else:
  r = 1.

for i in range(20):
  
  if r > -a:
    r = -a

  print("v=%f a=%f r=%r" % (v, a, r))

  v += a + r * 0.5
  a += r
 
 
