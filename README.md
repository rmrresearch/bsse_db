# Base DB

We can naively calculate the interaction energy between two monomers $A$ and $B$ using the following formula:

```math
\Delta E_{AB} = E_{AB} - E_A - E_B
```

However, this approach overestimates the interaction energy because the calculation of $E(AB)_{AB}$ artificially stabilizes monomers $A$ and $B$ by including the basis functions of the other monomer. This effect is known as **Basis Set Superposition Error (BSSE)**.

To correct for this, we calculate the energies of $A$ and $B$ in the **combined basis** (i.e., including ghost orbitals from the other monomer), resulting in the BSSE-corrected interaction energy:

```math
\Delta E(AB)_{AB} = E(AB)_{AB} - E(AB)_{A} - E(AB)_{B}
```

In this notation:

* The **parentheses** indicate the basis set used.
* The **subscript** indicates the monomers the energy is computed for.

So, the original (uncorrected) interaction energy becomes:

```math
\Delta E_{AB} = E(AB)_{AB} - E(A)_{A} - E(B)_{B}
```

The **BSSE** is the difference between the uncorrected and corrected interaction energies:

```math
\text{BSSE} = \Delta E(AB)_{AB} - \Delta E_{AB} = \sum_{i} \theta(AB)_i
```

Here, $\theta(AB)_i$ is the BSSE that monomer $i$ experiences when computed in the full $AB$ basis set rather than its own.

## Extension to Larger Systems

For larger systems, interaction energy and BSSE are defined analogously. In a three-body system, the **BSSE-free total interaction energy** 
includes both two-body and three-body terms:

```math
\text{Total Interaction Energy} = \sum_{i < j} \Delta E(ij)_{ij} + \Delta E(ABC)_{ABC}
```

Where:

```math
\Delta E(ij)_{ij} = E(ij)_{ij} - E(ij)_{i} - E(ij)_{j}
```

```math
\Delta E(ABC)_{ABC} = E(ABC)_{ABC} - \sum_i E(ABC)_i - \sum_{i < j} \Delta E(ABC)_{ij}
```

Note:
```math
\Delta E(ij)_{ij} \neq \Delta E(ABC)_{ij}
```
The former includes stabilization from the third basis set (i.e., it's not BSSE-free).

We define:

* $\Delta \epsilon_{IJ} = \Delta E_{IJ}(IJ)$ — the BSSE-free interaction energy
* $\Delta \theta_{IJ}(ABC) = \Delta E_{IJ}(ABC) - \Delta E_{IJ}(IJ)$ — the two-body BSSE from third-body stabilization

Using this, the total BSSE-free interaction energy simplifies to:

```math
\text{Total Interaction Energy} = E(ABC)_{ABC} - \sum_{i < j} \Delta \theta_{ij}(ABC) - \sum_i E(ABC)_i
```

## BSSE-Contaminated Interaction Energy

The interaction energy contaminated by BSSE is simply:

```math
\text{BSSE-Contaminated Interaction Energy} = E(ABC)_{ABC} - \sum_i E(i)_{i}
```

Taking the difference yields the **total BSSE**:

```math
\text{BSSE} = \sum_{i < j} \Delta \theta_{ij}(ABC) + \sum_i \theta_{i}(ABC)
```

Where:

```math
\theta_i(ABC) = E(i)_{ABC} - E(i)_{i}
```

So, the total BSSE consists of:

* **One-body BSSEs** $\theta_i$
* **Two-body BSSEs** $\Delta \theta_{ij}$

## Approximate BSSE Schemes

Computing all BSSE terms exactly requires evaluating every monomer in **every possible basis set** containing that monomer, which is computationally infeasible for large systems. This project evaluates various approximations for the BSSE.

1. **No monomer BSSE**:

   ```math
   \sum_i \theta_i(ABC) \approx 0
   ```

2. **Neglect third-body stabilization**:

   ```math
   \sum_i \theta_i(ABC) \approx \sum_i \theta_i(ij)
   ```

3. **No dimer BSSE**:

   ```math
   \sum_{i < j} \Delta \theta_{ij}(ABC) \approx 0
   ```
