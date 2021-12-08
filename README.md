# Li-ion-cell-capacity-estimation
Robust capacity estimation algorithm that can be incorporated into the BMSs of Li-ion battery systems

The battery capacity is essential for the applications such as cell balancing and SoC estimation.\
\
As the cell ages, the total capacity of the cell degrades. The lithium ion cells have side reactions that occur during charge and discharge cycles. These reactions consume lithium that could have been otherwise used during charge and discharge cycles, and structural deterioration of the electrode active materials that eliminates lithium storage sites.
The degradation over aging is called capacity fading and the objective of this project is to explore the algorithms that track capacity fade.
Approximate Weighted Total Least Squares estimation (based on [1]) is used and A test data set is taken from [2]\
\
•[1] Plett, G.L.,"Recursive Approximate Weighted Total Least Squares Estimation of Battery Cell Total Capacity" Journal of Power Sources 196(4), 2011 https://doi.org/10.1016/j.jpowsour.2010.09.048 \
•[2] Devie, Arnaud, George Baure, and Matthieu Dubarry. 2018. "Intrinsic Variability in the Degradation of a Batch of Commercial 18650 Lithium-Ion Cells" Energies 11, no. 5: 1031. 
https://doi.org/10.3390/en11051031 
