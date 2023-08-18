.numToString <- function(num, dec, p.val = FALSE) {
  #Function to convert numbers to strings with
  #specified number of decimal places
  string <- formatC(num, digits = dec, format = 'f')
  if (p.val & string < 0.001) {
    string <- '<0.001'
  }
  return(string)
}

.missingVar <- function(variables, df) {
  #Check if there are missing variables in
  #the list of variables of the baseline table
  missing.vars <- c()
  for (var in variables) {
    if (any(colnames(df) == var) == FALSE) {
      missing.vars <- append(missing.vars, var)
    }
  }
  return(missing.vars)
}

#Function to create a dataframe for factors
.varIsFactor <- function(factor, factor.label, perc.format, group,
                         group.levels, group.levels.n, df) {

  factor.levels <- levels(df[, factor])
  tbl.n <- table(df[, factor], df[, group])
  tbl.perc <- prop.table(tbl.n, perc.format)
  p.val <- .numToString(chisq.test(tbl.n)$p.value, 3, p.val = TRUE)

  #Empty rows for the output dataframe
  df.row <- data.frame()
  df.rows <- data.frame()

  for (i.factor in 1:length(factor.levels)) {

    for (i.group in 1:group.levels.n) {
      #Calculate numbers and percentages
      n <- tbl.n[factor.levels[i.factor], group.levels[i.group]]
      perc <- tbl.perc[factor.levels[i.factor], group.levels[i.group]] * 100
      #Create a formatted string
      n.perc.string <- paste0(.numToString(n, 0),
                              ' (',
                              .numToString(perc, 1),
                              ')')
      #Create a new dataframe from the formated string
      new.col <- data.frame(n.perc.string)
      #Column and row names
      colnames(new.col) <- group.levels[i.group]
      if (length(factor.levels) > 2) {
        rownames(new.col) <- paste0(factor.levels[i.factor])
      } else {
        rownames(new.col) <- paste0(factor.label,
                                    ' (',
                                    factor.levels[i.factor],
                                    ')')
      }
      #Check if its the first or last iteration in the loop
      if (i.group == 1) {
        df.row <- new.col
      } else if (i.group == length(group.levels)) {
        #Only add P value if the factor has 2 or less levels, otherwise p value
        #should only be in the heading
        p.val.col <- if (length(factor.levels) > 2) {data.frame(NA)} else {data.frame(p.val)}
        #Column and row names for p value
        colnames(p.val.col) <- ("P value")
        if (length(factor.levels) > 2) {
          rownames(p.val.col) <- paste0(factor.levels[i.factor])
        } else {
          rownames(p.val.col) <- paste0(factor.label,
                                        ' (',
                                        factor.levels[i.factor],
                                        ')')
        }
        df.row <- cbind(df.row, new.col, p.val.col)
      } else {
        df.row <- cbind(df.row, new.col)
      }
    }

    #Add to output dataframe
    df.rows = if(i.factor == 1) {

      #If there are more than 2 levels to a factor a heading is added to the
      #final dataframe with the name of the factor and a p value. In other cases
      #2 rows are created which include the rowname factor (level) and both of
      #the rows will contain a p value
      if (length(factor.levels) > 2) {
        for (i.group in 1:(group.levels.n + 1)) {
          if (i.group == 1) {
            header.row <- data.frame(NA)
            colnames(header.row)[1] <- group.levels[1]
          } else if (i.group == (group.levels.n + 1)) {
            p.val.col <- data.frame(p.val)
            colnames(p.val.col)[1] <- 'P value'
            header.row <- cbind(header.row, p.val.col)
            rownames(header.row)[1] <- factor.label
          } else {
            new.header.row <- data.frame(NA)
            colnames(new.header.row)[1] <- group.levels[i.group]
            header.row <- cbind(header.row, new.header.row)
          }
        }
      }

      if (length(factor.levels) > 2) {
        rbind(header.row,df.row)
      } else {
        df.row
      }
    } else {
      rbind(df.rows, df.row)
    }
  }
  return(df.rows)
}

.varIsNumeric <- function(numeric.var, numeric.label, distribution,
                          group, group.levels, group.levels.n, df) {
  #Calculate P value
  if (distribution == 'normal') {
    p.val <- if(group.levels.n == 2) {
      .numToString(
        t.test(df[, numeric.var] ~ df[, group])$p.value,
        3,
        p.val = TRUE)} else {
          .numToString(
            summary(aov(df[, numeric.var] ~ df[, group]))[[1]][5]$'Pr(>F)'[1],
            3,
            p.val = TRUE)
        }

    summ.mean <- aggregate(
      x = df[, numeric.var],
      by = list(df[, group]),
      FUN = mean,
      na.rm = TRUE
    )
    summ.sd <- aggregate(
      x = df[, numeric.var],
      by = list(df[, group]),
      FUN = sd,
      na.rm = TRUE
    )

  } else {
    p.val <- if(length(levels(df[, group])) == 2) {
      .numToString(wilcox.test(df[, numeric.var] ~ df[, group])$p.value,
                   3,
                   p.val = TRUE)} else {
                     .numToString(kruskal.test(df[, numeric.var] ~ df[, group])$p.value,
                                  3,
                                  p.val = TRUE)}

    summ.qnt <- aggregate(
      x = df[, numeric.var],
      by = list(df[, group]),
      FUN = quantile,
      na.rm = TRUE
    )

  }

  df.row <- data.frame()

  for (i.group in 1:group.levels.n) {

    level.name <- group.levels[i.group]

    if (distribution == "normal") {
      mean <- .numToString(summ.mean[which(summ.mean$Group.1 == level.name), ]$x,
                           1
      )
      sd <- .numToString(summ.sd[which(summ.sd$Group.1 == level.name), ]$x,
                         1
      )
      mean.sd <- paste0(mean, ' ± ', sd)
      new.col <- data.frame(mean.sd)

    } else {
      median <- .numToString(summ.qnt[which(summ.qnt$Group.1 == level.name), ]$x[, '50%'],
                             1
      )
      lower <- .numToString(summ.qnt[which(summ.qnt$Group.1 == level.name), ]$x[, '25%'],
                            1
      )
      upper <- .numToString(summ.qnt[which(summ.qnt$Group.1 == level.name), ]$x[, '75%'],
                            1
      )
      median.u.l <- paste0(median, ' [', lower, ' to ', upper, ']')
      new.col <- data.frame(median.u.l)
    }

    colnames(new.col) <- level.name
    rownames(new.col) <- paste0(numeric.label)

    if (i.group == 1) {
      df.row <- new.col
    } else if (i.group == group.levels.n) {
      p.val.col <- data.frame(p.val)
      colnames(p.val.col) <- ("P value")
      rownames(p.val.col) <- paste0(numeric.label)
      df.row <- cbind(df.row, new.col, p.val.col)
    } else {
      df.row <- cbind(df.row, new.col)
    }
  }
  return(df.row)
}


DescribeData.2 <- function(variables = variables,
                           normal = normal,
                           group = group,
                           df = df,
                           names = c(),
                           nmiss = TRUE,
                           perc.format = 2) {

  if (!is.null(.missingVar(variables, df))) {
    stop(paste0('Variables are missing from dataframe: ',
                .missingVar(variables, df)))
  }

  #Group levels and number of groups from the input group
  group.levels <- levels(df[, group])
  group.levels.n <- length(levels(df[, group]))

  n.variables <- length(variables)
  output.df <- data.frame()

  #Loop through all variables in the variables input string
  for (i.var in 1:n.variables) {
    var.name <- variables[i.var]

    var.label <- ifelse(length(names) != 0 & length(names) == n.variables,
                        names[i.var],
                        var.name)
    #Type of variable
    if (is.factor(df[, var.name])) {
      new.rows <- .varIsFactor(var.name, var.label, perc.format, group,
                               group.levels, group.levels.n, df)
    } else if (var.name %in% normal) {
      new.rows <- .varIsNumeric(var.name, var.label, 'normal', group,
                                group.levels, group.levels.n, df)
    } else {
      new.rows <- .varIsNumeric(var.name, var.label, 'skewed', group,
                                group.levels, group.levels.n, df)
    }

    if(i.var == 1) {
      output.df = new.rows
    } else {
      output.df <- rbind(output.df, new.rows)
    }

  }
  return(output.df)
}
