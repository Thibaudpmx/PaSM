lines <- readLines("D:/these/Second_project/QSP/modeling_work/VT_simeoni/calciumhomeostasis-boneresorption-model-master/model_c.txt")

# remove line of only comments
lines <- lines[!grepl("^ */", lines)]
# remove line of define
lines <- lines[!grepl("^ *#define ", lines)]
# comments line of section (starting with $)
lines <- gsub("^ *\\$",  "#$", lines)

# remove double
lines <- gsub("^ *double", "", lines)

# handle init sections

blocinistarts <- grep("\\$INIT", lines)
blociniends <- grep("\\$ODE ", lines)

new_bloc <- lines[blocinistarts:blociniends] %>%
  gsub(pattern = " *=",replacement =  "(0) =")

lines[blocinistarts:blociniends]  <- new_bloc

# remove ;

lines <- gsub(";", "", lines)

# handle power by removing pow and replaince "," by "**"
lines <- gsub("pow", "", lines)
lines <- gsub(",", "**", lines)

# new synthax diff eq

lineswithdxdt <- grep("dxdt_", lines)

new <- lines[lineswithdxdt] %>%
  gsub(pattern = "dxdt_", replacement = "d/dt(") %>%
  gsub(pattern = " *=", replacement = ") = ")

lines[lineswithdxdt] <- new

# comment lines with IPRED

lines[grep("CMTFLAG==", lines)] <- paste0("#", lines[grep("CMTFLAG==", lines)] )

# remove capture bloc
lines <- lines[1:(grep("\\$CAPTURE", lines) - 1)]

# remove // as comments
lines <- gsub("//", "#", lines)
lines <- gsub("/\\*","#", lines )
#add line jump
lines <- lines %>% paste0(collapse = "\n")

# Modif manual
#1) remove line useess (just a name of a cmt / parameter)
#2) line pkgut pkcent pkper1 pkper2 splitted
#3) all param with value declaration first, then param defined in relation to other param (just copy/past)

# add Rxode
parse_expr(paste0("RxODE({", lines, "})"))
eval()
paste0("R")
