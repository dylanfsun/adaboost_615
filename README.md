# adaboost_615
r package for biostat 615 adaboost project

Instructions:

1. install.packages("devtools")

2. library(devtools) 

3. install_github(repo = "dylanfsun/adaboost_615")

4. library(adaboost615)

5. parameters <- adaboost(test[1:400,])

6. predicted <- predictAda615(test[401:500, 1:1000], parameters)

7. error_rate(predicted, test[401:500, 1001])
