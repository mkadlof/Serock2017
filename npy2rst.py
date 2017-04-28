#!/usr/bin/env python3

import numpy as np

 def save_restraints(contact_map, domain):
 #    14 :35 red 0.2 300000.0
     contact_map = np.tril(contact_map)
     w = []
     X, Y = contact_map.nonzero()
     if len(X) != 0:
         for x,y in zip(X,Y):
             w.append(":{} :{} red 0.2 300000.0\n".format(x+1, y+1))
         w[-1] = w[-1][:-1]  # odciecie ostatniego przejscia do nowego wiersza
     else:
         with open('errors.txt','a') as err:
             err.write('brak restraints w ' + str(domain) + '\n')

     with open(os.path.join(domain.path, 'restraints.txt'), 'w') as restraints:
         restraints.writelines(w)