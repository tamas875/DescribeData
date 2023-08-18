\name{DescribeData}
\alias{DescribeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DescribeData
}
\description{
Make a baseline characteristics table.
}
\usage{
DescribeData(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{variables}{
List or vector of variables. Each variable should be given as a string.
}
\item{normal}{
List or vector of normally distributed variables. Each variable should be given as a string.
}
\item{group}{
Grouping variable (columns) as a string.
}
\item{names}{
List or vector of variable names. Each name should be given as a string.
}
\item{perc.format}{
How to format percentages for categorical variables. For more information see the prop.table documentation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
DescribeData (variables = variables, df = df)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
