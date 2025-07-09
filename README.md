# base_db
We can naively calculate the interaction energy between two monomers $A$ and $B$ with the following formula. 

```math
 \Delta E_{AB} = E_{AB} - E_A - E_B
```
\
This approach however overestimates the interaction energy because the calculation of $E(AB)\_{AB}$ artificial stabilizes the monomers of A with the inclusion of the B basis set (and likewise for B). This is refered to as Basis Set Suppperposition Error (BSSE).
\
To fix this we calculate the energy of both $A$ and $B$ in the combined basis to get the BSSE corrected interaction energy. 

```math
 \Delta E_{\text{corrected}}(AB)_{AB} = E(AB)_{AB} - E(AB)_{A} - E(AB)_{B}
```

The monomers in the parentheses refer to the basis set and the subscript refers to the monomers the energy is being calculated for. So in this notation the first equation, the interaction energy without BSSE correction, becomes

```math
 \Delta E(AB)_{AB} = E(AB)_{AB} - E(A)_A - E(B)_B
```
The BSSE is simply the difference between these two terms.
```math
\text{BSSE} = \Delta E(AB)_{AB} - \Delta E_{\text{corrected}}(AB)_{AB}
```
For larger systems the BSSE and Interaction Energy is defined analogously.
```math
 \Delta E_{\text{corrected}}(ABC)_{ABC} = E(ABC)_{ABC} - \Sigma_{i} E(ABC)_{i} + \Sigma_{i,j} (E(ij)_{ij} + E(ABC)_{ij})
```
```math
 \Delta E_(ABC)_{ABC} = E(ABC)_{ABC} - \Sigma_{i} E(i)_{i}
```
