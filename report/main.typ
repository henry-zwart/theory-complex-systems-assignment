#import "@local/uva-game-theory:0.1.0": assignment

#set text(lang: "ENG", region: "GB")

#let authors = (
  (name: "John F. Nash", id: "31415926"),
  (name: "Julia Robinson", id: "27182818"),
)

#let course = (
  lecturer: "Dr. Ulle Endriss",
  name: "Game Theory",
  code: "5314GATH6Y",
)

#let collaboration = [
  Replace this text with your collaboration paragraph. 
  Replace this text with your collaboration paragraph. 
  Replace this text with your collaboration paragraph. 
  Replace this text with your collaboration paragraph. 
  Replace this text with your collaboration paragraph. 
  Replace this text with your collaboration paragraph. 
]

#show: assignment.with(
  title: "Homework #1",
  course: course,
  collaboration: collaboration,
  authors: authors,
  font: "New Computer Modern",
)

// Display main content (solutions)
#include "solutions.typ"
