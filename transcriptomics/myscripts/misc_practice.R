
# Playing in R to understand basic R functions and objects####



x <- 5
students <- data.frame( #makes data frame
  name = c("A", "B", "C"), #c() combines everything in ()
  height = c(62,68, 72)
)

head(students) #shows preview
class(students) #shows its a data frame
str(students) #types of data, chr = character num = number

students$name[1] #read me the first entry of "name: column
students[1,2] #look at (column,row)