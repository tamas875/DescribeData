\name{DescribeData}
\alias{DescribeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DescribeData
}
\description{
Used to create a baseline characteristics table. 
}
\usage{
DescribeData(variables, normal, group, df, names, nmiss = TRUE, perc.format = 2, between.iqr = ' to ')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variables}{A list of variables within the dataframe to be considered for the baseline table. Each variable should be given as a string.}
  \item{normal}{A list of variables that are normally distributed and exist in the variables parameter.}
  \item{group}{The name of the grouping variable (column) used for stratification. Should be given as a string.}
  \item{df}{The dataframe where the variables are found. Variable names should correspond to column names in the dataframe.}
  \item{names}{A list of row names to be specified in the baseline table. List order should match the order of variables in the "variables" parameter. Each name should be a string.}
  \item{nmiss}{A logical value (TRUE or FALSE) specifying whether to include the number of missing values in the column/variable. (Default is TRUE)}
  \item{perc.format}{Formatting for percentages in categorical variables. Refer to the prop.table documentation for more information.}
  \item{between.iqr}{Optional. Specifies the character that is placed between the interquartile ranges. By default is set to ' to '.}
}
\details{
Create a dataframe with column names defined by the "group" variable and row names specified by the "names" variable. 
For normally distributed variables, calculate the mean plus the standard deviation (mean + sd). 
For non-normally distributed variables, specify the median plus the interquartile range (median + iqr). 
For categorical variables, specify the count (n) and the percentage of individuals in this group.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
Tamas Szili-Torok
}

\examples{
# Example 1: Basic Usage
# Create a baseline characteristics table for a dataset.
# Consider variables "age," "income," and "education" stratified by the "gender" column.
# Specify custom row names "Age", "Income" and "Education"
baseline_data <- DescribeData(
  variables = c("age", "income", "education"),
  normal = c("age"),
  group = "gender",
  df = your_data_frame,
  names = c("Age (years)", "Income ($)", "Education"),
  nmiss = FALSE
)
View(baseline_data)

\begin{table}[]
\begin{tabular}{lllll}
\hline
            & Male                      & Female                  & P value          &  \\ \hline
Age (years) & 43.0 ± 12.5               & 40.0 ± 12.5             & \textless{}0.001 &  \\
Income (\$) & 2005.0 {[}1805 to 2293{]} & 2045 {[}1745 to 2354{]} & 0.004            &  \\
Education   & 300 (70)                  & 400 (60)                & \textless{}0.001 &  \\ \hline
\end{tabular}
\end{table}

# Example 2: Customized Interquartile Range Separator
# Create a baseline characteristics table with a custom separator for interquartile ranges.
# Use variables "height," "weight," and "blood_pressure," stratified by "smoking_status."
# Specify row names and customize the interquartile range separator as " - ".
baseline_data2 <- DescribeData(
  variables = c("height", "weight", "blood_pressure"),
  normal = c("blood_pressure"),
  group = "smoking_status",
  df = your_data_frame,
  names = c("Height (m)", "Weight (kg)", "Systolic Blood Pressure (mmHg)"),
  between.iqr = ' - '
)
View(baseline_data2)

\begin{table}[]
\begin{tabular}{lllll}
\hline
                               & Never smoked        & Currently Smokes    & Formerly smoked     & P value          \\ \hline
Height (cm)                    & 176 {[}161 - 189{]} & 163 {[}158 - 176{]} & 174 {[}154 - 192{]} & \textless{}0.001 \\
Weight (kg)                    & 87 {[} 69 - 103{]}  & 65 {[}54 - 77{]}    & 86 {[}56 - 106{]}   & 0.004            \\
Systolic Blood Pressure (mmHg) & 123.3 ± 10.4        & 132 ± 11.0          & 135 ± 14.3          & \textless{}0.001 \\ \hline
\end{tabular}
\end{table}
}
