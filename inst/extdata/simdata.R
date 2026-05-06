## simdata.R
## Generate a synthetic dataset with the same flavour as the NHANES II
## example used in the ginteff Stata help, but well-conditioned so that
## logit/ologit fits converge cleanly on the side-by-side comparisons.
## Writes simdata.csv (read by both examples_ginteff.R and
## examples_ginteff.do).

set.seed(20260430)

n <- 2000

## Demographics
female <- rbinom(n, size = 1, prob = 0.51)
race   <- sample(1:3, n, replace = TRUE, prob = c(0.70, 0.15, 0.15))
                                    # 1 = white, 2 = black, 3 = other

age    <- round(runif(n, 25, 75), 1)
height <- round(rnorm(n, mean = 175 - 12*female, sd = 8), 1)
weight <- round(rnorm(n, mean = 70 + 0.30*(age - 50) + 0.45*(height - 170),
                      sd = 11), 1)

## ---- Binary "poor health" outcome -----------------------------------------
##  True DGP includes (female x race) interaction and a small
##  (age x height x weight) 3-way interaction.
zage <- (age    - 50) / 10
zhei <- (height - 170) / 10
zwei <- (weight - 70)  / 10

eta_h <- -1.0 +
         0.40 * zage + 0.10 * zwei - 0.05 * zhei +
         0.30 * female +
         0.50 * (race == 2) - 0.10 * (race == 3) +
         0.40 * female * (race == 2) - 0.20 * female * (race == 3) +
         0.06 * zage * zhei * zwei

p_h        <- plogis(eta_h)
health_2l  <- rbinom(n, 1, p_h)

## ---- Ordered 5-level health ---------------------------------------------
eta_o <- 0.50 * zage + 0.20 * female + 0.20 * female * zage
y_lat <- eta_o + rlogis(n)
health <- as.integer(cut(y_lat,
                         breaks = c(-Inf, -1.2, -0.4, 0.4, 1.2, Inf),
                         labels = 1:5))
                                  # 1 = excellent ... 5 = poor

## ---- Diabetes (used in some Stata examples) -----------------------------
p_d      <- plogis(-2.0 + 0.30 * zage - 0.20 * female + 0.10 * zage * female)
diabetes <- rbinom(n, 1, p_d)

dat <- data.frame(
    id        = seq_len(n),
    health_2l = health_2l,
    health    = health,
    diabetes  = diabetes,
    female    = female,
    race      = race,
    age       = age,
    height    = height,
    weight    = weight
)

write.csv(dat, "simdata.csv", row.names = FALSE)
cat("Wrote simdata.csv with ",
    nrow(dat), " obs and ", ncol(dat), " variables.\n", sep = "")
cat("Health distribution: \n")
print(table(dat$health))
cat("Race x female cells:\n")
print(table(dat$race, dat$female))
