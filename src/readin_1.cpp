/**
 * @file readin.cpp
 * @brief Generic Input Parsing Utility.
 *
 * This module provides a robust mechanism for reading typed data from an input
 * stream. It handles:
 * - Skipping whitespace and comments (lines starting with '#').
 * - Type safety and conversion (Integer, Double, String, Char, Boolean).
 * - Error reporting with line numbers.
 * - Debugging verbosity levels.
 *
 * @author Sathwik Bharadwaj, Andrei Ilyashenko, Arya Retheeshan (Refactoring)
 * @date 2025-11-24
 */

#include "petsc.h"
#include <cstring> // Needed for strcpy
#include <fstream> // Needed for file io
#include <iostream>
#include <sstream> // Needed for string conversions
using namespace std;

/**
 * @brief Ignores characters in the stream until a delimiter is reached.
 *
 * This function consumes characters from the input stream until the specified
 * delimiter character `c` is found. The delimiter itself is also consumed.
 * Unlike `istream::ignore`, this function does not require a maximum count,
 * making it suitable for skipping variable-length content.
 *
 * @param file Input file stream.
 * @param c    Delimiter character to stop at.
 */
void ignore_until(ifstream &file, char c) {
  char ch;
  while (file.get(ch) && ch != c)
    ;
}

/**
 * @brief Determines the current line number in the input file.
 *
 * This function calculates the line number corresponding to the current
 * "get" position of the file stream. It does this by rewinding to the
 * beginning and counting newlines up to the current position.
 *
 * @note This operation can be expensive for large files as it re-reads from the
 * start. It is primarily used for error reporting.
 *
 * @param file Input file stream.
 * @return int The 1-based line number.
 */
int determine_current_line_number(ifstream &file) {
  int pos = file.tellg();
  file.seekg(0, ios::beg);

  int linenum = 0;
  while (file.tellg() <= pos) {
    linenum++;
    ignore_until(file, '\n');
  }

  file.clear();
  file.seekg(pos, ios::beg);
  return linenum;
}

/**
 * @brief Reads a typed value from an input stream with error checking.
 *
 * This is the core utility for parsing simulation input files. It reads a
 * string token from the file, strips comments, and attempts to convert it to
 * the requested type.
 *
 * Supported Types:
 * - 'I' / 'i': Integer
 * - 'D' / 'd': Double precision floating point
 * - 'S' / 's': String (C-style char array)
 * - 'C' / 'c': Single character
 * - 'B' / 'b': Boolean (accepts numeric 0/1 or "true"/"false")
 *
 * Debug Levels (ndebug):
 * - 0: Silent (only errors).
 * - 1: Verbose (print read values and warnings).
 * - 2: Strict (treat warnings as errors).
 *
 * @param file        Input file stream.
 * @param data        Pointer to the destination variable (void*).
 * @param description Description of the parameter (for logging/errors).
 * @param type        Character code indicating the expected data type.
 * @param units       String describing units (optional, for logging).
 * @param ndebug      Debug verbosity level.
 * @return PetscErrorCode 0 on success.
 * @throws int Throws 1 on fatal error (e.g., EOF, conversion failure).
 */
PetscErrorCode readin(ifstream &file, void *data, const char *description,
                      char type, const char *units, int ndebug) {
  stringstream s;
  string str;
  PetscErrorCode ierr;

  // -------------------------------------------------------------------------
  // 1. Read Raw Token
  // -------------------------------------------------------------------------
  // Read the next whitespace-delimited token.
  file >> skipws >> str;

  // -------------------------------------------------------------------------
  // 2. Handle Comments
  // -------------------------------------------------------------------------
  // If the token starts with '#', it's a comment. Ignore the rest of the line
  // and read the next token.
  while (str[0] == '#') {
    ignore_until(file, '\n');
    file >> skipws >> str;

    // Check for unexpected EOF
    if (file.eof()) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD,
                      "ERROR reached end of file while scanning for input %s.",
                      description);
      CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD, "The input type was %c", type);
      CHKERRQ(ierr);
      throw(1);
    }
  }

  if (file.eof()) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "ERROR reached end of file while scanning for input %s.",
                       description);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "The input type was %c", type);
    CHKERRQ(ierr);
    throw(1);
  }

  // -------------------------------------------------------------------------
  // 3. Strip Inline Comments
  // -------------------------------------------------------------------------
  // Check if a comment starts mid-token (e.g., "5#comment").
  // We want to keep "5" and push "#comment" back into the stream.
  int comment_location = str.length();
  for (int i = 0; i < str.length(); i++) {
    if (str[i] == '#') {
      comment_location = i;
      break;
    }
  }

  if (comment_location != str.length()) {
    // Push back characters from the end up to the comment start
    for (int i = str.length() - 1; i >= comment_location; i--) {
      file.putback(str[i]);
    }
    // Truncate the string to remove the comment part
    str.erase(str.begin() + comment_location, str.end());
  }

  // -------------------------------------------------------------------------
  // 4. Type Conversion
  // -------------------------------------------------------------------------
  s << str;
  switch (type) {
  // --- Integer ---
  case 'I':
  case 'i':
    int ansi;
    s >> ansi;
    if (s.fail()) {
      cout << "ERROR on line ";
      cout.fill('0');
      cout.width(4);
      cout << determine_current_line_number(file)
           << " expected integer for input " << description
           << ", but input was \"" << str << "\". Conversion failed." << endl;
      throw(1);
    }
    if (ndebug > 0) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD, "line %d : %s  %d\n",
                      determine_current_line_number(file), description, ansi);
      CHKERRQ(ierr);
    }
    // Check for partial conversion warnings
    if (!s.eof() && ndebug > 0) {
      cout << "----> expected integer input, but input was \"" << str
           << "\". Partial conversion occurred." << endl;
      if (ndebug > 1)
        throw(1);
    }
    (*(int *)data) = ansi;
    break;

  // --- Double ---
  case 'D':
  case 'd':
    double ansd;
    s >> ansd;
    if (s.fail()) {
      cout << "ERROR on line ";
      cout.fill('0');
      cout.width(4);
      cout << determine_current_line_number(file)
           << " expected double for input " << description
           << ", but input was \"" << str << "\". Conversion failed." << endl;
      throw(1);
    }
    if (ndebug > 0) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD, "line %d : %s  %f\n",
                      determine_current_line_number(file), description, ansd);
      CHKERRQ(ierr);
    }
    if (!s.eof() && ndebug > 0) {
      cout << "----> expected double input, but input was \"" << str
           << "\". Partial conversion occurred." << endl;
      if (ndebug > 1)
        throw(1);
    }
    (*(double *)data) = ansd;
    break;

  // --- String ---
  case 'S':
  case 's':
    if (ndebug > 0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD, "line %d : %s  %s\n",
                         determine_current_line_number(file), description,
                         str.c_str());
      CHKERRQ(ierr);
    }
    {
      // Safe string copy with truncation to PETSC_MAX_PATH_LEN
      char *dest = (char *)data;
      size_t maxlen = PETSC_MAX_PATH_LEN;
      std::strncpy(dest, str.c_str(), maxlen - 1);
      dest[maxlen - 1] = '\0';
    }
    break;

  // --- Character ---
  case 'C':
  case 'c':
    if (ndebug > 0) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD, "line %d : %s  %c\n",
                      determine_current_line_number(file), description, str[0]);
      CHKERRQ(ierr);
    }
    if (str.length() > 1 && ndebug > 0) {
      cout << "----> expected char input, but input was \"" << str
           << "\". Partial conversion occurred." << endl;
      if (ndebug > 1)
        throw(1);
    }
    (*(char *)data) = str[0];
    break;

  // --- Boolean ---
  case 'B':
  case 'b':
    bool ansb;
    s >> ansb;
    if (s.fail()) {
      s.clear();
      s >> boolalpha >> ansb; // Try reading "true"/"false"
    }
    if (s.fail()) {
      cout << "ERROR on line ";
      cout.fill('0');
      cout.width(4);
      cout << determine_current_line_number(file)
           << " expected boolean for input " << description
           << ", but input was \"" << str << "\". Conversion failed." << endl;
      throw(1);
    }
    if (ndebug > 0) {
      ierr =
          PetscPrintf(PETSC_COMM_WORLD, "line %d : %s is %d\n",
                      determine_current_line_number(file), description, ansb);
      CHKERRQ(ierr);
    }
    if (!s.eof() && ndebug > 0) {
      cout << "----> expected boolean input, but input was \"" << str
           << "\". Partial conversion occurred." << endl;
      if (ndebug > 1)
        throw(1);
    }
    (*(bool *)data) = ansb;
    break;

  default:
    cout << "The type \"" << type
         << "\" was not recognized while reading "
            "for input "
         << description << ". Allowed types: I, D, S, C, B." << endl;
    throw(1);
  }
  return 0;
}
