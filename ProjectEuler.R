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
    while (i < input.number){
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
max <- 0
for (i in 13:1000){
    current <- substr(mega, i - 12, i)
    current <- DigitProduct(current)
    if (current > max){
        max <- current
    }
}
max


##9 There exists exactly one Pythagorean triplet for which a + b + c = 1000. Find the product abc

vals <- 1:1000
squares <- vals*vals
triplets <- c()
for (i in 1:1000){
    for (j in 1:1000){
        if ((squares[i] + squares[j]) %in% squares){
            triplets <- c(triplets, vals[i], vals[j], sqrt(squares[i] + squares[j]))
        }
    }
}

triplets
hit <- c()
for (i in 1:(length(triplets) / 3)){
    vals <- triplets[c(3*i - 2, 3*i - 1, 3*i)]
    if (sum(vals) == 1000){
        hit <- 3*i
        break
    }
}

triplets[hit] * triplets[hit - 1] * triplets[hit -2]

## 10
## Why didn't this update??

## 11 What is the greatest product of 4 adjacent (horizontal, vertical, diagonol) numbers in the grid below
numbers <- matrix(c
(08, 02, 22, 97, 38, 15, 00, 40, 00, 75, 04, 05, 07, 78, 52, 12, 50, 77, 91, 08,
49, 49, 99, 40, 17, 81, 18, 57, 60, 87, 17, 40, 98, 43, 69, 48, 04, 56, 62, 00,
81, 49, 31, 73, 55, 79, 14, 29, 93, 71, 40, 67, 53, 88, 30, 03, 49, 13, 36, 65,
52, 70, 95, 23, 04, 60, 11, 42, 69, 24, 68, 56, 01, 32, 56, 71, 37, 02, 36, 91,
22, 31, 16, 71, 51, 67, 63, 89, 41, 92, 36, 54, 22, 40, 40, 28, 66, 33, 13, 80,
24, 47, 32, 60, 99, 03, 45, 02, 44, 75, 33, 53, 78, 36, 84, 20, 35, 17, 12, 50,
32, 98, 81, 28, 64, 23, 67, 10, 26, 38, 40, 67, 59, 54, 70, 66, 18, 38, 64, 70,
67, 26, 20, 68, 02, 62, 12, 20, 95, 63, 94, 39, 63, 08, 40, 91, 66, 49, 94, 21,
24, 55, 58, 05, 66, 73, 99, 26, 97, 17, 78, 78, 96, 83, 14, 88, 34, 89, 63, 72,
21, 36, 23, 09, 75, 00, 76, 44, 20, 45, 35, 14, 00, 61, 33, 97, 34, 31, 33, 95,
78, 17, 53, 28, 22, 75, 31, 67, 15, 94, 03, 80, 04, 62, 16, 14, 09, 53, 56, 92,
16, 39, 05, 42, 96, 35, 31, 47, 55, 58, 88, 24, 00, 17, 54, 24, 36, 29, 85, 57,
86, 56, 00, 48, 35, 71, 89, 07, 05, 44, 44, 37, 44, 60, 21, 58, 51, 54, 17, 58,
19, 80, 81, 68, 05, 94, 47, 69, 28, 73, 92, 13, 86, 52, 17, 77, 04, 89, 55, 40,
04, 52, 08, 83, 97, 35, 99, 16, 07, 97, 57, 32, 16, 26, 26, 79, 33, 27, 98, 66,
88, 36, 68, 87, 57, 62, 20, 72, 03, 46, 33, 67, 46, 55, 12, 32, 63, 93, 53, 69,
04, 42, 16, 73, 38, 25, 39, 11, 24, 94, 72, 18, 08, 46, 29, 32, 40, 62, 76, 36,
20, 69, 36, 41, 72, 30, 23, 88, 34, 62, 99, 69, 82, 67, 59, 85, 74, 04, 36, 16,
20, 73, 35, 29, 78, 31, 90, 01, 74, 31, 49, 71, 48, 86, 81, 16, 23, 57, 05, 54,
01, 70, 54, 71, 83, 51, 54, 69, 16, 92, 33, 48, 61, 43, 52, 01, 89, 19, 67, 48), 20, 20, byrow = T)

## goes through rows
max <- 0
for (i in 1:20){
    for (j in 1:17){
        temp <- numbers[i, j] * numbers[i, j + 1] * numbers[i, j + 2] * numbers[i, j + 3]
        if (temp > max){
            max <- temp
        }
    }
}

## goes through columns
for (i in 1:17){
    for (j in 1:20){
        temp <- numbers[i, j] * numbers[i + 1, j] * numbers[i + 2, j] * numbers[i + 3, j]
        if (temp > max){
            max <- temp
        }
    }
}

## goes through diagonal
for (i in 1:17){
    for (j in 1:17){
        temp <- numbers[i, j] * numbers[i + 1, j + 1] * numbers[i + 2, j + 2] * numbers[i + 3, j + 3]
        if (temp > max){
            max <- temp
        }
    }
}
## goes through opposite diagnonal
for (i in 1:17){
    for (j in 1:17){
        temp <- numbers[i +3, j] * numbers[i + 2, j + 1] * numbers[i + 1, j + 2] * numbers[i, j + 3]
        if (temp > max){
            max <- temp
        }
    }
}





##12 The sequence of triangle numbers is generated by adding the natural numbers. 
## What is the value of the first triangle number to have over five hundred divisors?

## this took about 10 minutes to run

NumberDivisors2 <- function(triangle){
    ## to count 1 and itself
    number <- 2
    ## only proceede if it's divisble by common numbers: if it is, add those common numbers and their reciprocols
    if (triangle %% 2 == 0 & triangle %% 3 == 0 & triangle %% 4 == 0 & triangle %% 5 == 0 & 
        triangle %% 6 == 0 & triangle %% 7 == 0){
        number <- number + 12
        for (i in 8:((triangle / 7) - 1)){
            if (triangle %% i == 0){
                number <- number + 1
            }
        }
    }else{
        ## if not, won't be first to 500
        
    }
    return(number)
}
found <- F
tri <- 28
i <- 7
while (found == F){
    num <- NumberDivisors2(tri)
    if (num > 500){
        break
    }else{
        i <- i + 1
        tri <- tri + i
    }
    
}



##13 what are the first 10 digits of the sum of the following 100 50 digit numbers?
big <- c(37107287533902102798797998220837590246510135740250,
46376937677490009712648124896970078050417018260538,
74324986199524741059474233309513058123726617309629,
91942213363574161572522430563301811072406154908250,
23067588207539346171171980310421047513778063246676,
89261670696623633820136378418383684178734361726757,
28112879812849979408065481931592621691275889832738,
44274228917432520321923589422876796487670272189318,
47451445736001306439091167216856844588711603153276,
70386486105843025439939619828917593665686757934951,
62176457141856560629502157223196586755079324193331,
64906352462741904929101432445813822663347944758178,
92575867718337217661963751590579239728245598838407,
58203565325359399008402633568948830189458628227828,
80181199384826282014278194139940567587151170094390,
35398664372827112653829987240784473053190104293586,
86515506006295864861532075273371959191420517255829,
71693888707715466499115593487603532921714970056938,
54370070576826684624621495650076471787294438377604,
53282654108756828443191190634694037855217779295145,
36123272525000296071075082563815656710885258350721,
45876576172410976447339110607218265236877223636045,
17423706905851860660448207621209813287860733969412,
81142660418086830619328460811191061556940512689692,
51934325451728388641918047049293215058642563049483,
62467221648435076201727918039944693004732956340691,
15732444386908125794514089057706229429197107928209,
55037687525678773091862540744969844508330393682126,
18336384825330154686196124348767681297534375946515,
80386287592878490201521685554828717201219257766954,
78182833757993103614740356856449095527097864797581,
16726320100436897842553539920931837441497806860984,
48403098129077791799088218795327364475675590848030,
87086987551392711854517078544161852424320693150332,
59959406895756536782107074926966537676326235447210,
69793950679652694742597709739166693763042633987085,
41052684708299085211399427365734116182760315001271,
65378607361501080857009149939512557028198746004375,
35829035317434717326932123578154982629742552737307,
94953759765105305946966067683156574377167401875275,
88902802571733229619176668713819931811048770190271,
25267680276078003013678680992525463401061632866526,
36270218540497705585629946580636237993140746255962,
24074486908231174977792365466257246923322810917141,
91430288197103288597806669760892938638285025333403,
34413065578016127815921815005561868836468420090470,
23053081172816430487623791969842487255036638784583,
11487696932154902810424020138335124462181441773470,
63783299490636259666498587618221225225512486764533,
67720186971698544312419572409913959008952310058822,
95548255300263520781532296796249481641953868218774,
76085327132285723110424803456124867697064507995236,
37774242535411291684276865538926205024910326572967,
23701913275725675285653248258265463092207058596522,
29798860272258331913126375147341994889534765745501,
18495701454879288984856827726077713721403798879715,
38298203783031473527721580348144513491373226651381,
34829543829199918180278916522431027392251122869539,
40957953066405232632538044100059654939159879593635,
29746152185502371307642255121183693803580388584903,
41698116222072977186158236678424689157993532961922,
62467957194401269043877107275048102390895523597457,
23189706772547915061505504953922979530901129967519,
86188088225875314529584099251203829009407770775672,
11306739708304724483816533873502340845647058077308,
82959174767140363198008187129011875491310547126581,
97623331044818386269515456334926366572897563400500,
42846280183517070527831839425882145521227251250327,
55121603546981200581762165212827652751691296897789,
32238195734329339946437501907836945765883352399886,
75506164965184775180738168837861091527357929701337,
62177842752192623401942399639168044983993173312731,
32924185707147349566916674687634660915035914677504,
99518671430235219628894890102423325116913619626622,
73267460800591547471830798392868535206946944540724,
76841822524674417161514036427982273348055556214818,
97142617910342598647204516893989422179826088076852,
87783646182799346313767754307809363333018982642090,
10848802521674670883215120185883543223812876952786,
71329612474782464538636993009049310363619763878039,
62184073572399794223406235393808339651327408011116,
66627891981488087797941876876144230030984490851411,
60661826293682836764744779239180335110989069790714,
85786944089552990653640447425576083659976645795096,
66024396409905389607120198219976047599490197230297,
64913982680032973156037120041377903785566085089252,
16730939319872750275468906903707539413042652315011,
94809377245048795150954100921645863754710598436791,
78639167021187492431995700641917969777599028300699,
15368713711936614952811305876380278410754449733078,
40789923115535562561142322423255033685442488917353,
44889911501440648020369068063960672322193204149535,
41503128880339536053299340368006977710650566631954,
81234880673210146739058568557934581403627822703280,
82616570773948327592232845941706525094512325230608,
22918802058777319719839450180888072429661980811197,
77158542502016545090413245809786882778948721859617,
72107838435069186155435662884062257473692284509516,
20849603980134001723930671666823555245252804609722,
53503534226472524250874054075591789781264330331690)
small <- big / 10^ 35
total <- sum(small)

##14
## collatze sequence is defined by the follow rule: 14 n → n/2 (n is even) n → 3n + 1 (n is odd). 
## What is the longest sequence starting under 1 million
original.sequence <- 1:1000000
## Any number under 500,000 can be discarded, because you can get the same sequence by starting
## from the number that is twice as large. 
sequence <- 500000:1000000

## any number over 500k that fits the pattern of an odd number under 500k x 3 +1 can be discarded. 
# This is every other multiple of 3 + 1
pattern2 <- seq(500002, 1000000, 3)
sequence <- sequence[!(sequence %in% pattern2)]


## any number that fits the pattern (seq * 3 + 1) / 2 can be discarded
pattern3 <- sequence[!(sequence %% 2 == 0)]
pattern3 <- pattern3 * 3 + 1
pattern3 <- pattern3 / 2
sequence <- sequence[!(sequence %in% pattern3)]

collatz <- function(starting.value, counter){
    if (starting.value == 1){
        return(counter + 1)
    }else if (starting.value %% 2 == 0){
        collatz(starting.value / 2, counter + 1)
    }else{
        collatz((starting.value*3) + 1, counter + 1)
    }
}
max <- 0
idx <- 0
list <- c(0)
for (i in sequence){
    val <- collatz(i, 0)
    list <- c(list, val)
    if (val > max){
        max <- val
        idx <- i
    }
}
idx

## 837799, 525 loops

##15 




## 17 If the numbers 1 to 5 are written out in words: one, two, three, four, five, then there are 3 + 3 + 5 + 4 + 4 = 19 
## letters used in total. If all the numbers from 1 to 1000 (one thousand) inclusive were written out in words, 
## how many letters would be used?

## letters per hundred

## hundred and == 10
hundreds.place <- 10

## 20-90 twenty thirty fourty fifty sixty seventy eighty ninety == 48
normal.tens <- 48

## 1-9 one two three four five six seven eight nine == 36
ones.place <- 36





##18: What is the longest path possible starting at the top of the triangle ?

## collapse the triangle from the bottom up, replacing each value with one below it which is larger

r1 <- 75
r2 <- c(95, 64)
r3 <- c(17, 47, 82)
r4 <- c(18, 35, 87, 10)
r5 <- c(20, 04, 82, 47, 65)
r6 <- c(19, 01, 23, 75, 03, 34)
r7 <- c(88, 02, 77, 73, 07, 63, 67)
r8 <- c(99, 65, 04, 28, 06, 16, 70, 92)
r9 <- c(41, 41, 26, 56, 83, 40, 80, 70, 33)
r10 <- c(41, 48, 72, 33, 47, 32, 37, 16, 94, 29)
r11 <- c(53, 71, 44, 65, 25, 43, 91, 52, 97, 51, 14)
r12 <- c(70, 11, 33, 28, 77, 73, 17, 78, 39, 68, 17, 57)
r13 <- c(91, 71, 52, 38, 17, 14, 91, 43, 58, 50, 27, 29, 48)
r14 <- c(63, 66, 04, 68, 89, 53, 67, 30, 73, 16, 69, 87, 40, 31)
r15 <- c(04, 62, 98, 27, 23, 09, 70, 98, 73, 93, 38, 53, 60, 04, 23)

rows <- list(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15)
## takes in a pyramid, compares bottom two rows, adding value of greater pair in bottom row to second from bottom
## returns modified pyramid
CollapsePyramid <- function(pyramid){
    end <- length(pyramid)
    if (end == 1){
        return(pyramid[[1]][1])
    }else{
    upper <- pyramid[[end - 1]]
    lower <- pyramid[[end]]
        for (i in 1:length(upper)){
            if (lower[i] > lower[i + 1]){
                upper[i] <- upper[i] + lower[i]
            }else{
                upper[i] <- upper[i] + lower[i + 1]
            }
        }
    pyramid[[length(pyramid) -1]] <-  upper
    CollapsePyramid(pyramid[-length(pyramid)])
    }
}

CollapsePyramid(list(c(3), c(7,4), c(2,4,6), c(8,5,9,3), c(1,1,2,1,1)))

pyramid <- list(c(3), c(7,4), c(2,4,6), c(8,5,9,3))

CollapsePyramid(rows) ## == 1074