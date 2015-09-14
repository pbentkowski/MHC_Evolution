/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   newtestclass.cpp
 * Author: piotr
 *
 * Created on Sep 14, 2015, Sep 14, 2015 3:48:45 PM
 */

#include "newtestclass.h"
#include "H2Pinteraction.h"


CPPUNIT_TEST_SUITE_REGISTRATION(newtestclass);

newtestclass::newtestclass() {
}

newtestclass::~newtestclass() {
}

void newtestclass::setUp() {
}

void newtestclass::tearDown() {
}

void newtestclass::testH2Pinteraction() {
    H2Pinteraction h2Pinteraction();
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

void newtestclass::testPresentGeneAny() {
    genestring hostgene;
    genestring pathogene;
    int simil_mesure;
    H2Pinteraction h2Pinteraction;
    bool result = h2Pinteraction.presentGeneAny(hostgene, pathogene, simil_mesure);
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

void newtestclass::testPresentGeneRow() {
    genestring hostgene;
    genestring pathogene;
    int simil_mesure;
    H2Pinteraction h2Pinteraction;
    bool result = h2Pinteraction.presentGeneRow(hostgene, pathogene, simil_mesure);
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

void newtestclass::testDoesInfectedHeteroBetter() {
    Host& host;
    Pathogen& patho;
    int simil_mesure;
    H2Pinteraction h2Pinteraction;
    h2Pinteraction.doesInfectedHeteroBetter(host, patho, simil_mesure);
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

void newtestclass::testDoesInfectedHomoBetter() {
    Host& host;
    Pathogen& patho;
    int simil_mesure;
    H2Pinteraction h2Pinteraction;
    h2Pinteraction.doesInfectedHomoBetter(host, patho, simil_mesure);
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

void newtestclass::testDoesInfectedAllToAll() {
    Host& host;
    Pathogen& patho;
    int simil_mesure;
    H2Pinteraction h2Pinteraction;
    h2Pinteraction.doesInfectedAllToAll(host, patho, simil_mesure);
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

void newtestclass::testDoesInfectedHeteroOnePerSpec() {
    Host& host;
    Pathogen& patho;
    int simil_mesure;
    H2Pinteraction h2Pinteraction;
    h2Pinteraction.doesInfectedHeteroOnePerSpec(host, patho, simil_mesure);
    if (true /*check result*/) {
        CPPUNIT_ASSERT(false);
    }
}

