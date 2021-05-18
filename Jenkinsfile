#!groovy

pipeline {
  agent {
    label 'docker-runner'
  }
  environment {
    KHALEESI_SLACK_TOKEN = credentials('KHALEESI_SLACK_TOKEN')
    JENKINS_SLACKBOT = credentials('JENKINS_SLACKBOT')
  }
  stages{
    stage('Pull PKPDsim') {
      steps {
        echo 'Pulling PKPDsim'
        sh """
        sudo rm -rf PKPDsim2
        git clone git@github.com:InsightRX/PKPDsim2.git
        """
      }
    }
    stage('Run docker container') {
      environment {
        AWS_ACCESS_KEY_ID = credentials('AWS_ACCESS_KEY_ID')
        AWS_SECRET_ACCESS_KEY = credentials('AWS_SECRET_ACCESS_KEY')
      }
      steps {
        echo "Running container"
        sh """
        \$(aws ecr get-login --no-include-email --region us-west-2 &> /dev/null)
        docker build -t pkpdmap .
        docker run -d --name ${BUILD_TAG} pkpdmap:latest
        """
      }
    }
    stage('Run R CMD check') {
      steps {
        echo 'Checking PKPDmap'
        sh """
        docker exec -i ${BUILD_TAG} Rscript -e "devtools::check()"
        """
      }
    }
  }
  post {
    always {
      sh "docker rm -f ${BUILD_TAG} &>/dev/null && echo 'Removed container'"
    }
    failure {
      sh "chmod +x slack_notification.sh"
      sh "./slack_notification.sh"
    }
  }
}
