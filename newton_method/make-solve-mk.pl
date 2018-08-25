#!/usr/bin/perl

@least_square = ("through_pass",
                 "weighted",
                 "normal_equation",
                 "weighted_normal_equation");

@algorithm = ("PartialPivLU",
              "FullPivLU",
              "HouseholderQR",
              "ColPivHouseholderQR",
              "FullPivHouseholderQR",
              "LLT",
              "LDLT",
              "JacobiSVD");

print "OBJS_SOLVE = \\\n";

foreach $algo (@algorithm)
{
    foreach $ls (@least_square)
    {
        print "\tsolve-" . $ls . "-" . $algo . ".o \\\n";
        print "\tsolve_fast-" . $ls . "-" . $algo . ".o \\\n";
    }
}

print "\n";
