/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   newtestclass.h
 * Author: piotr
 *
 * Created on Sep 14, 2015, Sep 14, 2015 3:48:44 PM
 */

#ifndef NEWTESTCLASS_H
#define NEWTESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

class newtestclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(newtestclass);

    CPPUNIT_TEST(testH2Pinteraction);
    CPPUNIT_TEST(testPresentGeneAny);
    CPPUNIT_TEST(testPresentGeneRow);
    CPPUNIT_TEST(testDoesInfectedHeteroBetter);
    CPPUNIT_TEST(testDoesInfectedHomoBetter);
    CPPUNIT_TEST(testDoesInfectedAllToAll);
    CPPUNIT_TEST(testDoesInfectedHeteroOnePerSpec);

    CPPUNIT_TEST_SUITE_END();

public:
    newtestclass();
    virtual ~newtestclass();
    void setUp();
    void tearDown();

private:
    void testH2Pinteraction();
    void testPresentGeneAny();
    void testPresentGeneRow();
    void testDoesInfectedHeteroBetter();
    void testDoesInfectedHomoBetter();
    void testDoesInfectedAllToAll();
    void testDoesInfectedHeteroOnePerSpec();

};

#endif /* NEWTESTCLASS_H */

