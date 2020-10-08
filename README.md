# bioinfo_assignment2
### Bioinformatics - Global Alignment using linear gap penalty

### Technologies
- Java 1.8
- Maven

### Build Project
    mvn package
- Make sure version of Java running is 1.8 and Maven is installed

### Execute Jar
    java -jar target/bioinfo_assignment2-0.1.0.jar [sequence file] [substitution matrix file] [gap penalty]
- sequence file - file containing DNA or Protein sequence to align
- substitution matrix file - file containing the substitution matrix to use (i.e. BLOSUM62 for protein sequences)
- gap penalty value to be used in alignment 

#### Human-mouse P53 example
    mvn package
    java -jar target bioinfo_assignment2-0.1.0.jar human_mouseP53.txt BLOSUM62.txt -6
