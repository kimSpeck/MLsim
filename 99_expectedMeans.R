weights <- c(2, 2, 2)

dummyGrid_2x2 = expand.grid(D5 = c(1, -1), 
                            D6 = c(1, -1))

apply(dummyGrid_2x2 * weights, 1, sum)

dummyGrid_2x2[,"D5xD6"] <- dummyGrid_2x2["D5"] * dummyGrid_2x2["D6"]

apply(dummyGrid_2x2 * weights, 1, sum)
apply(dummyGrid_2x2 * c(2, -2, 2), 1, sum)
apply(dummyGrid_2x2 * c(2, 2, -2), 1, sum)
apply(dummyGrid_2x2 * c(-2, -2, 2), 1, sum)

dummyGrid_2x2x2 = expand.grid(D5 = c(1, -1), 
                        D6 = c(1, -1), 
                        D7 = c(1, -1))

apply(dummyGrid_2x2x2 * weights, 1, sum)