# AdjustedMutualInformationCpp

Adjusted Mutual Information (AMI) from  

Information Theoretic Measures for Clusterings Comparison: Is a Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J., in Procs. the 26th International Conference on Machine Learning (ICML'09)

which was implemented in Matlab by N.X. Vinh is now implemented in Cpp with a Matlab wrapper.

It obtains very good speed up. In particular when the number of clusters is large compared to the number of records.

E.g. N = 100 and two clusterings with 50 clusters: 

  AMImatlab = -0.034101 in 0.041557 seconds
  
  AMIcpp = -0.034101 in 0.000896 seconds
  
  Speedup = 46.3806 
