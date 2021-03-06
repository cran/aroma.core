%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaMicroarrayDataSetTuple.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AromaMicroarrayDataSetTuple}
\docType{class}
\alias{AromaMicroarrayDataSetTuple}


\title{The AromaMicroarrayDataSetTuple class}

\description{
 Package:  aroma.core \cr
\bold{Class AromaMicroarrayDataSetTuple}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSetList}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{AromaMicroarrayDataSetTuple}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.core]{AromaUnitTotalCnBinarySetTuple}}\cr

public abstract static class \bold{AromaMicroarrayDataSetTuple}\cr
extends \link[R.filesets]{GenericDataFileSetList}\cr


}

\usage{
AromaMicroarrayDataSetTuple(..., .setClass="AromaMicroarrayDataSet")
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSetList}}.}
  \item{.setClass}{The name of the class of the input set.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as} \tab  -\cr
 \tab \code{as.AromaMicroarrayDataSetTuple} \tab  -\cr
 \tab \code{getChipTypes} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getSets} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{\link[aroma.core:nbrOfChipTypes.AromaMicroarrayDataSetTuple]{nbrOfChipTypes}} \tab Gets the number of chip types.\cr
}


 \bold{Methods inherited from GenericDataFileSetList}:\cr
as, as.GenericDataFileSetList, as.character, as.data.frame, as.list, assertDuplicates, clone, extract, getAsteriskTags, getDefaultFullName, getFileList, getFileListClass, getFullNames, getNames, getSet, getSets, getTags, indexOf, length, nbrOfFiles, nbrOfSets, setTags

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
