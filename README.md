# Solar System
The goal of this project is to overcome the limitations of the analytical solution of the N-body problem, which does not admit a general closed-form solution, by implementing a numerical method to solve the equations of motion.

The physical system considered includes the Sun and the eight planets of the Solar System, each influenced by the gravitational interaction with all other bodies.

The equations of motion for each body are derived by applying Newton's Law of Universal Gravitation within the framework of Newton’s Second Law of Dynamics, accounting for the mutual interactions among all bodies in the system.

This leads to a system of coupled differential equations, which is solved numerically using the classical fourth-order Runge-Kutta (RK4) method.

Below is a brief overview of the implemented method:

### RK4 derivation method

      x_i+1 = x_i + dt*sum_i(b_i * k_i)
      ki = F(t + c_i*dt, x + dt*sum_j(a_ij * k_j)

      with:   F = F(t, x)
              a_ij, b_i, c_i from Butcher Tableau

      Butcher Tableau for RK4
      | c1: 0                                                     |
      | c2: 0.5      a21: 0.5                                     |
      | c3: 0.5      a31: 0       a32: 0.5                        |
      | c4: 1        a41: 0       a42: 0      a43: 1              |
      |              b1: 1/6      b2: 1/3     b3: 1/3     b4: 1/6 |

      so:      
      k1 = f( t, x )
      k2 = f( t + 0.5*dt, x + dt*(0.5*k1) )
      k3 = f( t + 0.5*dt, x + dt*(0.5*k2) )
      k4 = f( t + dt, x + dt*(1*k3) )

      x_i+1 = x_i + dt*[ (1/6 * k1) + (1/3 * k2) + (1/3 * k3) + (1/6 * k4) ] 
      
  Note: in our case, the system is time-invariant:
      f = f(x), not f = f(t, x)

### 

          



