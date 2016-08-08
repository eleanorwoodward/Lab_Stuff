## Project Euler Coding Attempts

##1: what is the sum of all multiples of 3 or 5 below 1000?

multiples.3 <- seq(3, 1000, 3)
multiples.5 <- seq(5, 999, 5)
keep.5 <- c(T, T, F)
new.multiples.5 <- multiples.5[keep.5]
sum(multiples.3, new.multiples.5)


##2 what is the sum of the even fibonacci numbers

FibGenerator <- function(max_value){
    fibs <- c(1,2)
    while (fibs[length(fibs)] < max_value){
        fibs <- c(fibs, sum(tail(fibs, 2)))
    }
    return(fibs[1:(length(fibs) - 1)])
}

keep.fib <- c(F,T,F)

sum(FibGenerator(4000000)[keep.fib])


##3 What is the largest prime factor of 600851475143
target <- 600851475143

Factorize <- function(input.number){
    i <- 2
    while (i < input.number){
        if (input.number %% i == 0){
            return(c(i, Factorize(input.number /  i)))
        }else{
            i <- i + 1
        }
    }
    return(input.number)
}

Factorize(target)


#4  A palindromic number reads the same both ways. The largest palindrome made from the product of two 2-digit numbers is 9009 = 91 × 99.
#Find the largest palindrome made from the product of two 3-digit numbers.

max <- 999*999

Place <- function(integer, num){
    
    ## calculate number of digits in the number
    num.dig <- ceiling(log10(integer))
    
    ## calculate number of tens place to shift our number over
    dig.shift <- num.dig - num 
    
    ## divide by 10 until the digit we want to return is in ones place
    val <- floor(integer / 10 ^ dig.shift)
    
    ## divide by 10, subtract difference, 
    lower.val <- val / 10
    return(10* (lower.val - floor(lower.val)))
}

candidates <- c()
for (i in 100:999){
    candidates <- c(candidates, sapply(100:999, function(x,y){x*y}, i))
}


keepers <- c()
for (i in 1:length(candidates)){
    val <- candidates[i]
    digits <- ceiling(log10(val))
    pally <- F
    for (j in 1:floor(digits/2)){
        if (abs(Place(val, j) - Place(val, digits + 1 - j)) < .001){
            pally <- T
        }else{
            pally <- F
            break
        }
    }
    if (pally){
        keepers <- c(keepers, val)
    }
}

max(keepers)


##5 What is the smallest number that is evenly divisible by the numbers 1:20
primes <- 3*5*7*11*13*17*19*16*3


##6 find the difference between the sum of the squares of 1:100 and the square of the sum of 1:100

x <- 1:100
sum(x*x) - (sum(x)*sum(x))

#7 starting with 2 as the first prime number, what is the 10,001st prime number

IsPrime <- function(input.number){
    i <- 2
    while (i <= sqrt(input.number)){
        if (input.number %% i == 0){
            return(F)
        }else{
            i <- i + 1
        }
    }
    return(T)
}

i <- 2
primes <- 0
while(primes < 10001){
    if (IsPrime(i)){
        primes <- primes + 1
    }
    i <- i + 1
}

##8 find the 13 adjacent numbers with the greatest product
mega <- "7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450"
substr(mega, 1, 13)
DigitProduct <- function(string){
    front <- substr(string, 1,1)
    rest <- substr(string, 2, nchar(string))
    if (rest == ""){
        return(as.numeric(front))
    }else{
        return(as.numeric(front) * DigitProduct(rest))
    }
}


## find the sum of all primes below 2 million

GeneratePrimes <- function(max_val){
    i <- 12
    primes <- c(2,3,5,7)
    while (i < max_val){
        val <- IsPrime(i - 1)
        if (val){
            primes <- c(primes, i - 1)
        }
        val <- IsPrime(i + 1)
        if (val){
            primes <- c(primes, i + 1)
        }
        i <- i + 6
    }
    return(primes)
}

fourmil <- GeneratePrimes(2000000)
y <- sum(fourmil)

## what is the largest product of four adjacent (vertical, horizontal, diaganal) numbers in the matrix?

values <- c(08 ,02 ,22 ,97 ,38 ,15 ,00 ,40 ,00 ,75 ,04 ,05 ,07 ,78 ,52 ,12 ,50 ,77 ,91 ,08,
49 ,49 ,99 ,40 ,17 ,81 ,18 ,57 ,60 ,87 ,17 ,40 ,98 ,43 ,69 ,48 ,04 ,56 ,62 ,00,
81 ,49 ,31 ,73 ,55 ,79 ,14 ,29 ,93 ,71 ,40 ,67 ,53 ,88 ,30 ,03 ,49 ,13 ,36 ,65,
52 ,70 ,95 ,23 ,04 ,60 ,11 ,42 ,69 ,24 ,68 ,56 ,01 ,32 ,56 ,71 ,37 ,02 ,36 ,91,
22 ,31 ,16 ,71 ,51 ,67 ,63 ,89 ,41 ,92 ,36 ,54 ,22 ,40 ,40 ,28 ,66 ,33 ,13 ,80,
24 ,47 ,32 ,60 ,99 ,03 ,45 ,02 ,44 ,75 ,33 ,53 ,78 ,36 ,84 ,20 ,35 ,17 ,12 ,50,
32 ,98 ,81 ,28 ,64 ,23 ,67 ,10 ,26 ,38 ,40 ,67 ,59 ,54 ,70 ,66 ,18 ,38 ,64 ,70,
67 ,26 ,20 ,68 ,02 ,62 ,12 ,20 ,95 ,63 ,94 ,39 ,63 ,08 ,40 ,91 ,66 ,49 ,94 ,21,
24 ,55 ,58 ,05 ,66 ,73 ,99 ,26 ,97 ,17 ,78 ,78 ,96 ,83 ,14 ,88 ,34 ,89 ,63 ,72,
21 ,36 ,23 ,09 ,75 ,00 ,76 ,44 ,20 ,45 ,35 ,14 ,00 ,61 ,33 ,97 ,34 ,31 ,33 ,95,
78 ,17 ,53 ,28 ,22 ,75 ,31 ,67 ,15 ,94 ,03 ,80 ,04 ,62 ,16 ,14 ,09 ,53 ,56 ,92,
16 ,39 ,05 ,42 ,96 ,35 ,31 ,47 ,55 ,58 ,88 ,24 ,00 ,17 ,54 ,24 ,36 ,29 ,85 ,57,
86 ,56 ,00 ,48 ,35 ,71 ,89 ,07 ,05 ,44 ,44 ,37 ,44 ,60 ,21 ,58 ,51 ,54 ,17 ,58,
19 ,80 ,81 ,68 ,05 ,94 ,47 ,69 ,28 ,73 ,92 ,13 ,86 ,52 ,17 ,77 ,04 ,89 ,55 ,40,
04 ,52 ,08 ,83 ,97 ,35 ,99 ,16 ,07 ,97 ,57 ,32 ,16 ,26 ,26 ,79 ,33 ,27 ,98 ,66,
88 ,36 ,68 ,87 ,57 ,62 ,20 ,72 ,03 ,46 ,33 ,67 ,46 ,55 ,12 ,32 ,63 ,93 ,53 ,69,
04 ,42 ,16 ,73 ,38 ,25 ,39 ,11 ,24 ,94 ,72 ,18 ,08 ,46 ,29 ,32 ,40 ,62 ,76 ,36,
20 ,69 ,36 ,41 ,72 ,30 ,23 ,88 ,34 ,62 ,99 ,69 ,82 ,67 ,59 ,85 ,74 ,04 ,36 ,16,
20 ,73 ,35 ,29 ,78 ,31 ,90 ,01 ,74 ,31 ,49 ,71 ,48 ,86 ,81 ,16 ,23 ,57 ,05 ,54,
01 ,70 ,54 ,71 ,83 ,51 ,54 ,69 ,16 ,92 ,33 ,48 ,61 ,43 ,52 ,01 ,89 ,19 ,67 ,48)

matrix.v <- matrix(values, 20, 20, byrow = T)

## check horizontal

max.h <- 0
loc.h <- c()
for (h in 1:20){    
    for (i in 1:17){
        val <- (matrix.v[h, i:(i+3)])
        product <- val[1] * val[2] * val[3] * val[4]
        if (product > max.h){
            max.h <- product
            loc.h <- c(h,i)
        }
    }
}

## check vertical
max.v <- 0
loc.v <- c()
for (h in 1:20){    
    for (i in 1:17){
        val <- (matrix.v[i:(i+3), h])
        product <- val[1] * val[2] * val[3] * val[4]
        if (product > max.v){
            max.v <- product
            loc.v <- c(i, h)
        }
    }
}


## diag
max.d <- 0
loc.d <- c()
for (h in 1:17){    
    for (i in 1:17){
        val <- c(matrix.v[h, i], matrix.v[h + 1, i + 1], matrix.v[h + 2, i + 2], matrix.v[h + 3, i + 3])
        product <- val[1] * val[2] * val[3] * val[4]
        if (product > max.d){
            max.d <- product
            loc.d <- c(i, h)
        }
    }
}



