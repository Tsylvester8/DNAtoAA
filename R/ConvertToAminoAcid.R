#' Convert a given DNA sequence or a RNA sequence to an animo acid sequence
#' @param x a vector of character strings
#' @return  The input DNA or RNA sequence will be converted to an amino acid sequence
#' @examples
#' x <- "AATGTTAT"
#' ConvertToAminoAcid(x)
#' @export



ConvertToAminoAcid <- function(x){
  if(nchar(x) %% 3 != 0){
    warning("The DNA/RNA string is not a multiple of 3")
  }

  input <- toupper(x) # get the upper case of the input string

  input <- gsub("T", "U", input)


  codonTable <- matrix(data = NA,
                       nrow = 64,
                       ncol = 2)

  num.codons <- nchar(input) / 3


  pos.codon <- seq(from = 1, by = 3, length.out = num.codons)

  AAseq <- c()


  codonTable[,1] <- c("UUU",
                      "UUC",
                      "UUA",
                      "UUG",
                      "UCU",
                      "UCC",
                      "UCA",
                      "UCG",
                      "UAU",
                      "UAC",
                      "UAA",
                      "UAG",
                      "UGU",
                      "UGC",
                      "UGA",
                      "UGG",
                      "CUU",
                      "CUC",
                      "CUA",
                      "CUG",
                      "CCU",
                      "CCC",
                      "CCA",
                      "CCG",
                      "CAU",
                      "CAC",
                      "CAA",
                      "CAG",
                      "CGU",
                      "CGC",
                      "CGA",
                      "CGG",
                      "AUU",
                      "AUC",
                      "AUA",
                      "AUG",
                      "ACU",
                      "ACC",
                      "ACA",
                      "ACG",
                      "AAU",
                      "AAC",
                      "AAA",
                      "AAG",
                      "AGU",
                      "AGC",
                      "AGA",
                      "AGG",
                      "GUU",
                      "GUC",
                      "GUA",
                      "GUG",
                      "GCU",
                      "GCC",
                      "GCA",
                      "GCG",
                      "GAU",
                      "GAC",
                      "GAA",
                      "GAG",
                      "GGU",
                      "GGC",
                      "GGA",
                      "GGG")

  codonTable[,2] <- c("F",
                      "F",
                      "L",
                      "L",
                      "S",
                      "S",
                      "S",
                      "S",
                      "Y",
                      "Y",
                      "*",
                      "*",
                      "C",
                      "C",
                      "*",
                      "W",
                      "L",
                      "L",
                      "L",
                      "L",
                      "P",
                      "P",
                      "P",
                      "P",
                      "H",
                      "H",
                      "Q",
                      "Q",
                      "R",
                      "R",
                      "R",
                      "R",
                      "I",
                      "I",
                      "I",
                      "M",
                      "T",
                      "T",
                      "T",
                      "T",
                      "N",
                      "N",
                      "K",
                      "K",
                      "S",
                      "S",
                      "R",
                      "R",
                      "V",
                      "V",
                      "V",
                      "V",
                      "A",
                      "A",
                      "A",
                      "A",
                      "D",
                      "D",
                      "E",
                      "E",
                      "G",
                      "G",
                      "G",
                      "G")

  for(i in pos.codon){

    codon <- substr(input, i, i+2)

    if(codon %in% codonTable[,1]){
      hit <- which(codon == codonTable[,1])

      AAseq <- c(AAseq, codonTable[hit, 2])
    }
  }
  return(AAseq)
}




