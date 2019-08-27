#!groovy

  pipeline {
    agent {
      label 'r-slave'
    }
    stages{
      stage('Dependencies - build json2test') {
        steps {
          echo 'building json2test'
          sh """
            #!/bin/bash
            set -ex
            pwd
            sudo chmod 777 ~/workspace
            cd ~/workspace
            if [ -d "json2test" ]; then
              sudo rm -R json2test
            fi
            git clone git@github.com:InsightRX/json2test.git
            cd json2test
            chmod +x slack_notification.sh
            R CMD INSTALL . --library=/usr/lib/R/site-library || { export STATUS=failed
            ./slack_notification.sh
            exit 1
            }
          """
        }
      }
      stage('Dependencies - build clinPK') {
        steps {
          echo 'building clinPK'
          sh """
            cd /$workspace
            if [ -d "clinPK2" ]; then
              sudo rm -R clinPK2
            fi
            git clone git@github.com:InsightRX/clinPK2.git
            cd clinPK2
            chmod +x slack_notification.sh
            R CMD INSTALL . --library=/usr/lib/R/site-library || { export STATUS=failed
            ./slack_notification.sh
            exit 1
            }
            """
        }
      }
      stage('Dependencies - PKPDsim') {
        steps {
          echo 'building PKPDsim'
          sh """
          if [ -d "PKPDsim2" ]; then
            sudo rm -R PKPDsim2
          fi
          git clone git@github.com:InsightRX/PKPDsim2.git
          cd PKPDsim2
          chmod +x slack_notification.sh
          R CMD INSTALL . --library=/usr/lib/R/site-library || { export STATUS=failed
          ./slack_notification.sh
          exit 1
          }
          R CMD check . --no-manual || { export STATUS=failed
          ./slack_notification.sh
          exit 1
          }
          """
        }

      }
      stage('Build - PKPDmap') {
        steps {
          echo 'buildilng PKPDmap'
          sh """
          if [ -d "PKPDmap" ]; then
            sudo rm -R PKPDmap
          fi
          git clone git@github.com:InsightRX/PKPDmap.git
          cd PKPDmap
          git checkout $GIT_BRANCH
          git pull origin $GIT_BRANCH
          chmod +x slack_notification.sh
          { R CMD INSTALL . --library=/usr/lib/R/site-library

          R CMD check . --no-manual
          } || { export STATUS=failed
          ./slack_notification.sh
          exit 1
          }
          """
        }
      }
    }
  }