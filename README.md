# Sugar Removal Utility (SRU)
##### An algorithmic approach for <i>in silico</i> removal of circular and linear sugars from molecular structures

## Description
Here, we present source code and examples for the Sugar Removal Utility, an algorithmic approach for <i>in silico</i> 
removal of circular and linear sugars from molecular structures, as described in [**put publication reference**]. 
<br>The algorithm's implementation is available in three forms: As a web application, a command-line application, and
as source code readily usable for other software development projects. Every form is open and free to use. The web 
application is available at [https://sugar.naturalproducts.net](https://sugar.naturalproducts.net) and its source code 
can be found [here](https://github.com/mSorok/SugarRemovalWeb). The command-line application along with its source code 
and the sugar removal algorithm's main implementation are part of this repository.
<p>Further description on the implemented sugar removal algorithm and its various configurations will be added soon.

## Contents
### SugarRemovalUtility CMD App
The sub-folder "SugarRemovalUtility CMD App" contains the sugar removal command-line application downloadable as 
compressed archive. After decompression, the JAR file "SugarRemovalUtility-jar-with-dependencies.jar" can be executed
from the command-line using Java version 11 or higher. A detailed explanation how to use the application can be found in
"Usage instructions.txt". Also, an example input file is provided, named "smiles_test_file.txt".

### Natural product test set
The text file "hand_picked_np.txt" contains a list of SMILES codes serving as a natural product test set for the 
performance of the SugarRemovalUtility. They were hand-picked from public databases via the 
[COlleCtion of Open NatUral producTs (COCONUT)](https://coconut.naturalproducts.net). More details can be found in the test class (see below) and the Sugar Removal 
Utility publication [**put reference**].

### Sources
In the directory <i>/src/main/java/de/unijena/cheminf/deglycosylation/</i> the class <i>SugarRemovalUtility</i> can be found.
This class represents the stand-alone implementation of the sugar removal algorithm. It can be used to detect and remove 
circular and linear sugar moieties from molecules supplied as CDK IAtomContainer objects with many configurable options.
Further documentation can be found in its JavaDoc comments.
The other sources available in <i>/src/main/java/de/unijena/cheminf/deglycosylation/</i> belong to the command-line 
application.

The class <i>SugarRemovalUtilityTest</i> can be found in the directory 
<i>/src/test/java/de/unijena/cheminf/deglycosylation/</i>. It is a JUnit test class that tests the performance of the 
Sugar Removal Utility on multiple specific molecular structures of natural products hand-picked from public databases 
(see above). Code examples of how to use and configure the <i>SugarRemovalUtility</i> class can be found here.

## Installation
### Command-line application JAR
The command-line application JAR has to be downloaded and decompressed. After that, it can be executed from the command-line
as described in the usage instructions. Java version 11 or higher has to be installed on your machine.

### Sources
This is a Maven project. In order to use the source code for your own software, download or clone the repository and 
open it in a Maven-supporting IDE (e.g. IntelliJ) as a Maven project and execute the pom.xml file. Maven will then take
care of installing all dependencies.

## Dependencies
* Java Development Kit (JDK) version 11
    * [AdoptOpenJDK](https://adoptopenjdk.net) (as one possible source of the JDK)
* Chemistry Development Kit (CDK) version 2.3
    * [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* Apache Maven version 4
    * [Apache Maven](http://maven.apache.org)
* JUnit version 4.13
    * [JUnit 4](https://junit.org/junit4/)

## References and useful links
**Sugar Removal Utility**
* [**put publication reference**]
* [Sugar Removal Web Application](https://sugar.naturalproducts.net)
* [Source Code of Web Application](https://github.com/mSorok/SugarRemovalWeb)

**Chemistry Development Kit (CDK)**
* [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* [Steinbeck C, Han Y, Kuhn S, Horlacher O, Luttmann E, Willighagen EL. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J Chem Inform Comput Sci. 2003;43(2):493-500.](https://dx.doi.org/10.1021%2Fci025584y)
* [Steinbeck C, Hoppe C, Kuhn S, Floris M, Guha R, Willighagen EL. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr Pharm Des. 2006; 12(17):2111-2120.](https://doi.org/10.2174/138161206777585274)
* [May JW and Steinbeck C. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014; 6:3.](https://dx.doi.org/10.1186%2F1758-2946-6-3)
* [Willighagen EL, Mayfield JW, Alvarsson J, Berg A, Carlsson L, Jeliazkova N, Kuhn S, Pluska T, Rojas-Chertó M, Spjuth O, Torrance G, Evelo CT, Guha R, Steinbeck C, The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J Cheminform. 2017; 9:33.](https://doi.org/10.1186/s13321-017-0220-4)
* [Groovy Cheminformatics with the Chemistry Development Kit](https://github.com/egonw/cdkbook)

**COlleCtion of Open NatUral producTs (COCONUT)**
* [COCONUT home page](https://coconut.naturalproducts.net)
* [Sorokina, M., Steinbeck, C. Review on natural products databases: where to find data in 2020. J Cheminform 12, 20 (2020).](https://doi.org/10.1186/s13321-020-00424-9)
