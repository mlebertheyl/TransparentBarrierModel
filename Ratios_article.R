#we´re using
prior.range[,1]
#corr for point 9
corr.dng[[90]][[1]][[9]]
#coords
id.coord9 <- matrix(c(id.node.dng[[49]]$id.coord[[9]][1], 
                      id.node.dng[[49]]$id.coord[[9]][2]), 
                      ncol = 2)

df9 <- df4plot(corr.dng[[90]][[1]][[9]], 
               dims = 300,
               id.coord = id.coord9)


####################################################################
#what's inside of function corr.plots (Ratios.Rmd)
df <- df9 %>% dplyr::filter(field >= 0.2)
#ggplot fiel vs dist
ggplot(data = df, aes(x = field, y = dist)) +
  geom_line()
#ggplot dist vs field
ggplot(data = df, aes(x = dist, y = field)) +
  geom_line()
#df summary
summary(df)
min.field <- df %>% dplyr::filter(field == min(df$field))
#what's the correlation at 21 dist, prior range 21
dist21 <- df9 %>% 
  dplyr::filter(dist >= 21) %>% 
  dplyr::filter(dist <= 22)
range21 <- dist21 %>% dplyr::filter(dist == min(dist21$dist))
#then range is the distance at which correlation is approx 0.13
#this is the matrix with the distance at which correlation is ~0.13 in the barriers depending on the range fraction scenario
x.sel21 <- range21$dist*x.sel

#for my "new definition" of range which I will call transparency
df0.8 <- df9 %>% dplyr::filter(field >= 0.8)
summary(df0.8)
min.field0.8 <- df0.8 %>% dplyr::filter(field == min(df0.8$field))
#the new range is the distance at which corr is 0.8
df0.9 <- df9 %>% dplyr::filter(field >= 0.9)
min.field0.9 <- df0.9 %>% dplyr::filter(field == min(df0.9$field))


####################################################################
#find the derivative of the curve
df0.13 <- df9 %>% dplyr::filter(field >= 0.13)
ggplot(data = df0.13, aes(x = field, y = dist)) +
  geom_line()
writexl::write_xlsx(df0.13, 'df0.13.xlsx')

spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
dy <- predict(spline_fit, deriv = 1)  # first derivative

par(mfrow=c(2,2))
plot(df0.13$field, df0.13$dist, main = "Data and Derivative")
plot(dy$x, dy$y, col = "blue")  # derivative line

# Find the derivative at x = ...
predict(spline_fit, x = 0.13, deriv = 1)$y
predict(spline_fit, x = 0.8, deriv = 1)$y
#curve for the corr is more or less linear at 0.5:
x0 <- 0.5
derivative_at_0.5 <- predict(spline_fit, x = x0, deriv = 1)$y
#let's say 

#1.Directly Modifying the Derivative
#If you want the derivative values themselves to decrease 70% faster than the original, you could scale the derivative by a factor.
#For "decay 70% faster," typically this means the new slope at any x is 70% larger in magnitude (steeper drop) than before.

x0 <- 0.5
original_deriv <- predict(spline_fit, x = x0, deriv = 1)$y
# If you want the derivative to decay 70% faster (i.e., 1.7 times steeper):
faster_decay_deriv <- original_deriv * 1.7
print(faster_decay_deriv)

#2.Adjusting the Spline's Smoothing or the Function
#If you want the whole curve to decay faster, you’d need to change your underlying model. For a smoothing spline, you could:
#Increase the penalty for wiggliness (increase spar in smooth.spline)
#Transform your x values to make the curve "compress" horizontally.
#If x goes to x / 1.7, everything happens 1.7 times "faster" along the x-axis.

spline_fit_fast <- smooth.spline(df0.13$field / 1.7, df0.13$dist)
deriv_fast <- predict(spline_fit_fast, x = x0 / 1.7, deriv = 1)$y
print(deriv_fast)

#Which method you should use depends on your modeling intent—whether you want a faster-decaying slope everywhere, or a faster-changing curve overall
#option 2 even though the result is the same for 0.5
plot(spline_fit_fast$x, spline_fit_fast$y, col = "blue")

# Fit smoothing spline to x compressed by 1.7 (70% faster decay)
#spline_fit_fast <- smooth.spline(df0.13$field / 1.7, df0.13$dist)

# Find the y value at x = 0.13 (original scale, so input 0.13/1.7 to the new model)
x_new <- 0.13
y_at_0.13_fast <- predict(spline_fit_fast, x = x_new)$y
print(y_at_0.13_fast)

x0 <- 0.13
y_at_0.13 <- predict(spline_fit, x = x0)$y
print(y_at_0.13)

#this is the ratio at 0.13
y_at_0.13_fast/y_at_0.13

#Find the y value at x = 0.2 (original scale, so input 0.2/1.7 to the new model)
x_new <- 0.2
y_at_0.2 <- predict(spline_fit_fast, x = x_new / 1.7)$y
print(y_at_0.2)
x0 <- 0.2
y_at_0.2_fast <- predict(spline_fit_fast, x = x0)$y
print(y_at_0.2_fast)
y_at_0.2_fast/y_at_0.2

#Find the y value at x = 0.5 (original scale, so input 0.5/1.7 to the new model)
x_new <- 0.5
y_at_0.5 <- predict(spline_fit_fast, x = x_new / 1.7)$y
print(y_at_0.5)
x0 <- 0.5
y_at_0.5_fast <- predict(spline_fit_fast, x = x0)$y
print(y_at_0.5_fast)
y_at_0.5_fast/y_at_0.5

#I need the ratio y_at_0.13_fast/y_at_0.13 to be 0.2
#scaling of smooth.spline(df0.13$field / 1.7, df0.13$dist) 
#such that y_at_0.13_fast/y_at_0.13 gives me 0.2

#find the scaling factor s so that
#predict(smooth.spline(df0.13$field / s, df0.13$dist), x = 0.13)$y / predict(spline_fit, x = 0.13)$y == 0.2.
#This can be solved by writing a small function and using uniroot() to search for s that achieves this ratio.

#Define a function that returns the difference between the ratio and 0.2
ratio_fun <- function(s) {
  spline_fit_fast <- smooth.spline(df0.13$field / s, df0.13$dist)
  y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
  ratio <- y_at_0.13_fast / y_at_0.13
  return(ratio - 0.2)
}

# Use uniroot to solve for s such that ratio == 0.2
# (choose an interval that makes sense for your case, e.g. 1.01 to 20)
result <- uniroot(ratio_fun, interval = c(1.01, 20))
s_solution <- result$root
cat("Scaling factor s needed:", s_solution, "\n")

# Check the result
spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
final_ratio <- y_at_0.13_fast / y_at_0.13

cat("Final ratio (should be 0.2):", final_ratio, "\n")

# Check the result
spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
final_ratio <- y_at_0.13_fast / y_at_0.13

######################################################################
#Define a more general function
#Generalized ratio function
ratio_fun <- function(s, x, aim) {
  spline_fit_fast <- smooth.spline(df0.13$field / s, df0.13$dist)
  y_at_x_fast <- predict(spline_fit_fast, x = x)$y
  y_at_x <- predict(spline_fit, x = x)$y
  ratio <- y_at_x_fast / y_at_x
  return(ratio - aim)
}

#Define target x and aim
target_x <- 0.5
aim_ratio <- 0.2

#Find s such that the ratio at target_x is aim_ratio
result <- uniroot(function(s) ratio_fun(s, target_x, aim_ratio), interval = c(1.01, 20))
s_solution <- result$root
cat("Scaling factor s needed:", s_solution, "\n")

# Check the result
spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
y_at_x_fast <- predict(spline_fit_fast, x = target_x)$y
y_at_x <- predict(spline_fit, x = target_x)$y
final_ratio <- y_at_x_fast / y_at_x

cat("Final ratio (should be", aim_ratio, "):", final_ratio, "\n")

#####################################################################
#loc.corr <- matrix(c(-909, 2875, -930.4101, 2900.680, -935.7784, 2899.498), ncol = 2, byrow = TRUE)
id.coord <- list()
for (l in c(2,8,9)) {
  id.coord[[l]] <- matrix(c(id.node.dng[[49]]$id.coord[[l]][1], 
                            id.node.dng[[49]]$id.coord[[l]][2]), 
                          ncol = 2)
}

df2 <- df4plot(corr.dng[[90]][[1]][[2]], 
               dims = 300,
               id.coord = id.coord[[2]])
df2 <- df2 %>% dplyr::filter(field >= 0.13)
ggplot(data = df2, aes(x = field, y = dist)) +
  geom_line()

#generalized function so I can pass df
ratio_fun <- function(df, s, x, aim) {
  # Fit original spline
  spline_fit <- smooth.spline(df$field, df$dist)
  
  # Fit scaled spline
  spline_fit_fast <- smooth.spline(df$field / s, df$dist)
  
  # Get predictions at x
  y_at_x <- predict(spline_fit, x = x)$y
  y_at_x_fast <- predict(spline_fit_fast, x = x)$y
  
  # Compute ratio
  ratio <- y_at_x_fast / y_at_x
  return(ratio - aim)
}

#Define target x and aim
target_x <- 0.13
aim_ratio <- 0.1
#Find s such that the ratio at target_x is aim_ratio
result2 <- uniroot(function(s) ratio_fun(df2, s, target_x, aim_ratio), interval = c(1.01, 20))
s_solution2 <- result2$root
cat("Scaling factor s needed:", s_solution2, "\n")

# Check the result
# (Optional, not needed in the function but useful for confirmation)
spline_fit2 <- smooth.spline(df2$field, df2$dist)
spline_fit_fast2 <- smooth.spline(df2$field / s_solution2, df2$dist)
y_at_x2 <- predict(spline_fit2, x = target_x)$y
y_at_x_fast2 <- predict(spline_fit_fast2, x = target_x)$y
final_ratio2 <- y_at_x_fast2 / y_at_x2
cat("Final ratio (should be", aim_ratio, "):", final_ratio2, "\n")

################################################################################
#
s_solution2 <- list()
x.tar <- c(0.13, 0.5, 0.8)
a.sol <- seq(0.1, 0.9, 0.1)
#point 2
for (x in 1:length(x.tar)) {
  s_solution2[[x]] <- list()
  for (a in 1:length(a.sol)) {
    target_x <- x.tar[x]
    aim_ratio <- a.sol[a]
    result2 <- uniroot(function(s) ratio_fun(df2, s, target_x, aim_ratio), interval = c(1.01, 20))
    s_solution2[[x]][[a]] <- result2$root
  }
}
#point 9
s_solution9 <- list()
for (x in 1:length(x.tar)) {
  s_solution9[[x]] <- list()
  for (a in 1:length(a.sol)) {
    target_x <- x.tar[x]
    aim_ratio <- a.sol[a]
    result2 <- uniroot(function(s) ratio_fun(df0.13, s, target_x, aim_ratio), interval = c(1.01, 20))
    s_solution9[[x]][[a]] <- result2$root
  }
}

#point 8
df8 <- df4plot(corr.dng[[90]][[1]][[8]], 
               dims = 300,
               id.coord = id.coord[[8]])
df8 <- df8 %>% dplyr::filter(field >= 0.13)
ggplot(data = df8, aes(x = field, y = dist)) +
  geom_line()

s_solution8 <- list()
for (x in 1:length(x.tar)) {
  s_solution8[[x]] <- list()
  for (a in 1:length(a.sol)) {
    target_x <- x.tar[x]
    aim_ratio <- a.sol[a]
    result2 <- uniroot(function(s) ratio_fun(df0.13, s, target_x, aim_ratio), interval = c(1.01, 20))
    s_solution8[[x]][[a]] <- result2$root
  }
}
#data frame
d <- a.sol
#point 2
a <- unlist(s_solution9[[1]])
s <- unlist(s_solution9[[2]]) #and 1-b
#point 9
aa <- unlist(s_solution2[[1]])
ss <- unlist(s_solution2[[2]])
dd <- unlist(s_solution2[[3]])
#point 8
aaa <- unlist(s_solution8[[1]])
sss <- unlist(s_solution8[[2]])
ddd <- unlist(s_solution8[[3]])

scales.df <- data.frame(d, a, s, aa, ss)#, aaa, sss)
data.frame(aa,ss,dd)
# Fit original spline and scaled spline (if not already available)
spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
spline_fit_fast <- smooth.spline(df0.13$field / 1.865426, df0.13$dist)

# Get derivatives at x_new
x_new <- 0.13
deriv_original <- predict(spline_fit, x = x_new, deriv = 1)$y
deriv_scaled   <- predict(spline_fit_fast, x = x_new, deriv = 1)$y

# Calculate the ratio (scaled/original)
deriv_ratio <- deriv_scaled / deriv_original

cat("Ratio of scaled to original derivative at x =", x_new, "is", deriv_ratio, "\n")

deriv_ratio_fun <- function(df, s, x) {
  spline_fit <- smooth.spline(df$field, df$dist)
  spline_fit_fast <- smooth.spline(df$field / s, df$dist)
  
  deriv_original <- predict(spline_fit, x = x, deriv = 1)$y
  deriv_scaled   <- predict(spline_fit_fast, x = x, deriv = 1)$y
  
  ratio <- deriv_scaled / deriv_original
  return(ratio)
}

deriv_ratio_fun(df0.13, 1.865426, 0.13)

# Define a sequence of x values for plotting (covering the observed range)
x_seq <- seq(min(df0.13$field), max(df0.13$field), length.out = 200)

# Compute derivatives for each x in the sequence
deriv_original <- predict(spline_fit, x = x_seq, deriv = 1)$y
deriv_scaled   <- predict(spline_fit_fast, x = x_seq, deriv = 1)$y

# Plotting
plot(x_seq, deriv_original, type = "l", lwd = 2, col = "blue",
     ylab = "Derivative of dist", xlab = "field",
     main = "Derivative Curves: Original vs. Scaled")
lines(x_seq, deriv_scaled, lwd = 2, col = "red")
legend("topright", legend = c("Original Derivative", "Scaled Derivative"),
       col = c("blue", "red"), lwd = 2)


# Plot the original data and spline_fit
plot(df0.13$field, df0.13$dist, 
     xlab = "field", ylab = "dist", 
     main = "Original Data and Scaled Data", 
     pch = 16, col = "blue", cex = 0.7)
# Add original spline fit
lines(spline_fit$x, spline_fit$y, col = "blue", lwd = 2)

# Plot the compressed data and spline_fit_fast
points(df0.13$field / s_solution, df0.13$dist, 
       pch = 17, col = "red", cex = 0.7)
# Add scaled spline fit
lines(spline_fit_fast$x, spline_fit_fast$y, col = "red", lwd = 2)

legend("topright", 
       legend = c("Original data + spline", "Scaled data + spline"),
       col = c("blue", "red"), pch = c(16,17), lwd = 2)

###############################################################################
# Fit original spline and scaled spline (if not already available)
spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
spline_fit_fast <- smooth.spline(df0.13$field / 6.058733, df0.13$dist)

# Get derivatives at x_new
x_new <- 0.5
deriv_original <- predict(spline_fit, x = x_new, deriv = 1)$y
deriv_scaled   <- predict(spline_fit_fast, x = x_new, deriv = 1)$y

# Calculate the ratio (scaled/original)
deriv_ratio <- deriv_scaled / deriv_original

cat("Ratio of scaled to original derivative at x =", x_new, "is", deriv_ratio, "\n")

deriv_ratio_fun(df0.13, 6.058733, x_new)

# Define a sequence of x values for plotting (covering the observed range)
x_seq <- seq(min(df0.13$field), max(df0.13$field), length.out = 200)

# Compute derivatives for each x in the sequence
deriv_original <- predict(spline_fit, x = x_seq, deriv = 1)$y
deriv_scaled   <- predict(spline_fit_fast, x = x_seq, deriv = 1)$y

# Plotting
plot(x_seq, deriv_original, type = "l", lwd = 2, col = "blue",
     ylab = "Derivative of dist", xlab = "field",
     main = "Derivative Curves: Original vs. Scaled")
lines(x_seq, deriv_scaled, lwd = 2, col = "red")
legend("topright", legend = c("Original Derivative", "Scaled Derivative"),
       col = c("blue", "red"), lwd = 2)


# Plot the original data and spline_fit
plot(df0.13$field, df0.13$dist, 
     xlab = "field", ylab = "dist", 
     main = "Original Data and Scaled Data", 
     pch = 16, col = "blue", cex = 0.7)
# Add original spline fit
lines(spline_fit$x, spline_fit$y, col = "blue", lwd = 2)

# Plot the compressed data and spline_fit_fast
points(df0.13$field / s_solution, df0.13$dist, 
       pch = 17, col = "red", cex = 0.7)
# Add scaled spline fit
lines(spline_fit_fast$x, spline_fit_fast$y, col = "red", lwd = 2)

legend("topright", 
       legend = c("Original data + spline", "Scaled data + spline"),
       col = c("blue", "red"), pch = c(16,17), lwd = 2)

###############################################################################
# Assuming you have:
s_solution <- 1.862555
# df0.13 already loaded

# Get the original and scaled spline fits
spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)

# Find predictions at x = 0.13
y_at_0.13 <- predict(spline_fit, x = 0.13)$y
y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y

# Compute the ratio
achieved_ratio <- y_at_0.13_fast / y_at_0.13
cat("Achieved ratio at x = 0.13 for s =", s_solution, "is", achieved_ratio, "\n")

for (s in 1:length(ss)) {
  # Assuming you have:
  s_solution <- ss[s]
  # df0.13 already loaded
  
  # Get the original and scaled spline fits
  spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
  spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
  
  # Find predictions at x = 0.13
  y_at_0.13 <- predict(spline_fit, x = 0.13)$y
  y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
  
  # Compute the ratio
  achieved_ratio <- y_at_0.13_fast / y_at_0.13
  print(cat("Achieved ratio at x = 0.13 for s =", s_solution, "is", achieved_ratio, "\n"))
}

for (s in 1:length(dd)) {
  # Assuming you have:
  s_solution <- dd[s]
  # df0.13 already loaded
  
  # Get the original and scaled spline fits
  spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
  spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
  
  # Find predictions at x = 0.13
  y_at_0.13 <- predict(spline_fit, x = 0.13)$y
  y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
  
  # Compute the ratio
  achieved_ratio <- y_at_0.13_fast / y_at_0.13
  print(cat("Achieved ratio at x = 0.13 for s =", s_solution, "is", achieved_ratio, "\n"))
}

for (s in 1:length(aa)) {
  # Assuming you have:
  s_solution <- aa[s]
  # df0.13 already loaded
  
  # Get the original and scaled spline fits
  spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
  spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
  
  # Find predictions at x = 0.13
  y_at_0.13 <- predict(spline_fit, x = 0.13)$y
  y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
  
  # Compute the ratio
  achieved_ratio <- y_at_0.13_fast / y_at_0.13
  print(cat("Achieved ratio at x = 0.13 for s =", s_solution, "is", achieved_ratio, "\n"))
}


for (s in 1:length(sss)) {
  # Assuming you have:
  s_solution <- sss[s]
  # df0.13 already loaded
  
  # Get the original and scaled spline fits
  spline_fit <- smooth.spline(df0.13$field, df0.13$dist)
  spline_fit_fast <- smooth.spline(df0.13$field / s_solution, df0.13$dist)
  
  # Find predictions at x = 0.13
  y_at_0.13 <- predict(spline_fit, x = 0.13)$y
  y_at_0.13_fast <- predict(spline_fit_fast, x = 0.13)$y
  
  # Compute the ratio
  achieved_ratio <- y_at_0.13_fast / y_at_0.13
  print(cat("Achieved ratio at x = 0.13 for s =", s_solution, "is", achieved_ratio, "\n"))
}
