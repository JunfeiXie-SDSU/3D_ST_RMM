# 3-D Smooth Turn Random Mobility Models for Aerial Vehicles
Matlab code for the two 3-D Smooth Turn (ST) random mobility models for aerial vehicles introduced in **J. Xie, Y. Wan, B. Wang, S. Fu, K. Lu, J. H. Kim, “A Comprehensive 3-D Modeling Framework for Airborne Networks”, IEEE Access, Vol. 6, pp. 22849-22862, 2018.**

The first model (z-independent ST) assumes that the z-dimensional movement is independent from the movement on the horizontal plane. This model captures normal aerial maneuvers, mostly observed in civilian applications. 

The second model (z-dependent ST) captures the correlation of movement among all three dimensions. This model can describe military aircraft performing 
"high-g" turns. 

Both models capture the correlation of accelerations across temporal and spatial directions during turns that comply with physical laws.

## Instruction 
To run the z-indepent ST model, run the `z_independent_withboundary.m` file. 

To run the z-dependent ST model, run the `z_dependent_withboundary.m` file.

## Paper citation
Please cite the following paper if you used the code or any of the models.

@article{xie2018comprehensive,
  title={A comprehensive 3-dimensional random mobility modeling framework for airborne networks},
  author={Xie, Junfei and Wan, Yan and Wang, Baoqian and Fu, Shengli and Lu, Kejie and Kim, Jae H},
  journal={IEEE Access},
  volume={6},
  pages={22849--22862},
  year={2018},
  publisher={IEEE}
}
```
