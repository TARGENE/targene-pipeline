# Understanding TarGene's Outputs

The successful execution of a TarGene run will result in the creation of a `results` folder containing at least 3 files:

- results.hdf5
- results.json
- QQ.png

## The QQ.png

The QQ.png is what you would expect, a Q-Q plot of your results. Since in TarGene you can use multiple estimators at once, this Q-Q plot will show all of them like the following.

![GWAS_QQ](../assets/gwas_QQ.png)

## The results.json

The `results.json` contains summary statistics of your run and will likely be all you need.

## The results.hdf5

The `results.hdf5` contains the extensive list of estimation results.